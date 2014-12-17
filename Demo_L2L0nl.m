
clear variables, close all, fclose('all'); clc

%% add paths ==============================================================

addpath(genpath('Data')); %HARDI data (fiber crossing)
addpath(genpath('NIFTI')); %nifti functions
addpath(genpath('Misc')); %miscellanious functions
addpath(genpath('Solvers')); %solvers
addpath(genpath('Operators')); %proximity operators
addpath(genpath('spgl1-1.7')); %SPGL1 library http://www.cs.ubc.ca/~mpf/spgl1/


%% define signal characteristics ==========================================

b			= 2000;
nSAMPLES	= 15;
SNR			= 10;
nDIRECTIONS = 200; % num. of rotations of the kernel (to generate the dictionary)

%% constants and flag definitions =========================================

data_path='Data/Phantom_Data/';
estimdic_path='Data/Phantom_Data/Estimated_Dicos';
maskfile= fullfile(data_path,'ISBI12__mask.nii');
fibersGTfile=fullfile(data_path,sprintf('ISBI12__dir.nii'));
nfibersGTfile=fullfile(data_path,'ISBI12__nfibers.nii');
fieldfile=fullfile(data_path,sprintf('ISBI12_%d_%ddB.mat',nSAMPLES,SNR));
bvecsfile=fullfile(data_path,sprintf('bvecs_%d.txt', nSAMPLES));
bvalsfile=fullfile(data_path,sprintf('bvals_%d.txt', nSAMPLES));
dirfile=fullfile(sprintf('gendir/gendir_%d.txt',nDIRECTIONS));
dicfile=fullfile(estimdic_path,sprintf('DICTIONARY_%d_%d_%ddB.mat',nDIRECTIONS, nSAMPLES, SNR));

SAVE_RESULTS=0;


%% Load files =============================================================

% mask -- voxels to recover the configuration from -- 
niiMASK = load_untouch_nii( maskfile );
    fprintf( '\t- MASK: %dx%dx%d with datatype %d (%d bits/voxel)\n', niiMASK.hdr.dime.dim(2), ...
        niiMASK.hdr.dime.dim(3), niiMASK.hdr.dime.dim(4), niiMASK.hdr.dime.datatype, niiMASK.hdr.dime.bitpix );
MASK=niiMASK.img; %mask
MASK=logical(MASK);


% ground truth fiber configuration (evaluation purposes) ----
% directions
niiGT  = load_untouch_nii( fibersGTfile );
    fprintf( '\t- PEAKS (GROUND TRUTH) : %dx%dx%dx%d with datatype %d (%d bits/voxel)\n', niiGT.hdr.dime.dim(2), ...
        niiGT.hdr.dime.dim(3), niiGT.hdr.dime.dim(4), niiGT.hdr.dime.dim(5), niiGT.hdr.dime.datatype, niiGT.hdr.dime.bitpix );
FIBERS_gt=niiGT.img;

% num of fibres per voxel
niiNFIBERS  = load_untouch_nii( nfibersGTfile );
    fprintf( '\t- nPEAKS (GROUND TRUTH) : %dx%dx%d with datatype %d (%d bits/voxel)\n', niiNFIBERS.hdr.dime.dim(2), ...
        niiNFIBERS.hdr.dime.dim(3), niiNFIBERS.hdr.dime.dim(4), niiNFIBERS.hdr.dime.datatype, niiNFIBERS.hdr.dime.bitpix );
nFIBERS_gt=niiNFIBERS.img;
%-------------------------------------------------------------

% field of directions and simulated signal:
% FIELD -> MultiTensor object with information of the configuration of the field of voxels
% SIGNAL -> 4-D tensor (nX, nY, nZ, nSAMPLES) simulated acquisition of nSAMPLES per voxel
%   with snr=SNR
load(fieldfile);

% sampling patern (bvecs and bvals) 
BVECS		= importdata( bvecsfile);
BVALS		= importdata( bvalsfile );


%% Files to store results ===================================================

filename=fullfile(sprintf('Demo_L2L0global_%ddB_%dsam.mat',SNR,nSAMPLES));
niifname=fullfile(sprintf('Demo_L2L0global_%ddB_%dsam.nii',SNR,nSAMPLES));
globalfilename=fullfile('L2L0global.mat');

%% Initialize variables: ==================================================

niiK=5; %number of expected directions in the niifile

nVOXELS=numel(MASK); % total number of voxels
actVOXELS=sum(MASK(:)); % total number of "active" voxels

nX	= size(FIELD,1);
nY	= size(FIELD,2);
nZ	= size(FIELD,3);

%% Simulate the signal ========================================================

XYZB=BVECS';XYZB=XYZB(2:end,:); % xyzb size (nSAMPLES,4)
XYZB(:,4) = BVALS(2:end);

%% Load dictionary ==========================================================

DIRECTIONS		= SphereGeneric( 'FromFile', dirfile, true );

load (dicfile);
nATOMS=size(DICTIONARY.ATOMs,2);

%% generate neighbourhood for each DIRECTION 

% NEIGH_DIR -> cell array of dimension (nDIRECTIONS,1) used to discriminate
%   the final fiber directions though a maxima search procedure
deg_max=15; 
[ NEIGH_DIR, DISTmin ] = neigh_dir_setup( DIRECTIONS, deg_max );

%% generate neighbourhood for each ATOM of the dictionary 

% NEIGH_ATOMS -> cell array of dimension (nATOMS,1) used to compute the weights

NEIGH_ATOMS=cell(nATOMS,1);

for d=1:nDIRECTIONS

    neigh=NEIGH_DIR{d};
    neigh=neigh(1:end/2);
    NEIGH_ATOMS{d}=[neigh;d];

end

d=nDIRECTIONS+1; %if the dictionary also contains isotropic compartments
while(d<=nATOMS)
    NEIGH_ATOMS{d}=[d];
    d=d+1;
end


%% generate neighbourhood for each VOXEL of the field

% vector (of dimension actVOXELS), for each active voxel it contains its index in the global nVOXELS 3D-vol.
allV_id=find(MASK(:)); 

actV_id=zeros(size(MASK));
for i=1:numel(MASK(:))
    if (find(allV_id==i)), actV_id(i)=find(allV_id==i); end
end

% generate neighborhood (radius r=1 voxel)
r=1;
NEIGH_VOX=neigh_vox_setup(ones(nX,nY,nZ),r);

% index translation
[ act_NEIGH_VOX ] = map_neigh_vox( NEIGH_VOX, allV_id, actV_id );


%% RECOVER the fiber configuration in each voxel
%  =============================================

% Initialize variables: =====================================================


% X \in R^{d \times av} (d: num directions//av: num active voxels)
X=zeros(nATOMS,actVOXELS); 

% nii structure to store fibers information (to be saved as nii file and
% read in software as Fiber Navigator)
niiFIBERS=niiMASK;
niiFIBERS.img=zeros(nX,nY,nZ, niiK*3); %all fiber directions per voxel 
niiFIBERS.hdr.dime.dim(1:5) = [4 nX nY nZ 3*niiK];
niiFIBERS.hdr.dime.bitpix = 32;
niiFIBERS.hdr.dime.datatype = 16;

FIBERS=zeros(nX,nY,nZ,niiK*3); % (niiK fiber directions per voxel)* 3 components/direction 
nFIBERS=zeros(nX,nY,nZ); %number of fibers per voxel

%resize and normalize measurements
IMG=reshape(SIGNAL, nX*nY*nZ, nSAMPLES+1);
IMG=bsxfun (@rdivide, double(IMG), double(IMG(:,1))); %normalize b0
IMG=IMG(:,2:end);
Y=IMG'; Y=Y(:,MASK(:));% define measurement matrix of dimension (nSAMPLES, actVOXELS)


M=actVOXELS*nSAMPLES;

%% Minimisation problem ==========


%% Define parameters

    %general parameters
    param.verbose=1;
    param.plot=0;
    param.rel_obj=1e-3;
    param.max_iter=1000;
    param.A=@(x)Phi(x,size(DICTIONARY.ATOMs,2),actVOXELS,DICTIONARY.ATOMs);
    param.At=@(x)Phit(x,size(DICTIONARY.ATOMs.',2),actVOXELS,DICTIONARY.ATOMs.');
    X_prev=zeros(nATOMS, actVOXELS);
    param.solini=X_prev(:);
    param.gamma=1/norm(DICTIONARY.ATOMs,2)^2;

    %L1 min parameters
    param.weights=1;
    param.Psi=@(x)x;param.Psit=@(x)x;param.nu_L1=1;
    param.pos_L1=0;
    param.tight_L1=1;
    param.max_iter_L1=1000;
    param.rel_obj_L1=1e-3;
    param.K=actVOXELS*3;

%% Solve minimisation problem

t0=tic;
    % 1st iteration
    [X, info]=solve_LASSO_FISTA(Y, param);
    X=reshape(X,nATOMS,actVOXELS);

    % smooth version for X
    [ Xs ] = fast_compute_weights(X, NEIGH_ATOMS, NEIGH_VOX, MASK);

    % fill the gaps when actVOXELS < nVOXELS
    X_complete=zeros(nATOMS,nVOXELS); % whole 3D-vol
    Xs_complete=zeros(nATOMS,nVOXELS);

    for i=1:actVOXELS

        X_complete(:,allV_id(i))=X(:,i);
        Xs_complete(:,allV_id(i))=Xs(:,i);

    end



%% Reweighting process 
    
    % parameters reweighting process
    tau=var(Xs(:));
    tau_thr=1e-7;
    param.plot=0;
    max_iter_rw=10;



    for iter_rw = 1:max_iter_rw

        %update parameters
        X_prev= X;
        W=1./(tau+abs(Xs));
        W(201,:)=1000; % no special weight for the isotropic compartment
        
        param.weights=W(:);
        param.solini=X(:);
        
        %solve minimisation problem
        [X, info]=solve_LASSO_FISTA(Y, param);

        X=reshape(X,nATOMS,actVOXELS);
        [ Xs ] = fast_compute_weights(X, NEIGH_ATOMS, NEIGH_VOX, MASK);

        rel_error=norm(X(:)-X_prev(:))/norm(X_prev(:));
        
        %check convergence
        if rel_error <= 1e-3
            break;
        end
        
        %update tau
        tau=max(tau/10,tau_thr);

    end
    
    
%% Compute final directions

    [ FIBERS_rw, nFIBERS_rw ] = compute_directions( X, niiMASK, NEIGH_DIR, DIRECTIONS, niiK );

    time=toc(t0);
    
%% Compute quality metrics

    [ SCORE ] = quality_measures( FIBERS_rw, FIBERS_gt,nFIBERS_gt, niiMASK);
    
%% Print results

    fprintf( '\n');
    fprintf( '   ===== L2L0_nl: Evaluation of the reconstruction ====================== \n');
    fprintf( '\n');
    fprintf( '\t- DATA VOLUME: %dx%dx%d voxels\n', niiMASK.hdr.dime.dim(2), niiMASK.hdr.dime.dim(3), niiMASK.hdr.dime.dim(4));
    fprintf( '\t- Mean Angular Error (over all voxels): %.2f\n', mean(SCORE.Ea(:)));
    fprintf( '\t- Median Angular Error (over all voxels): %.2f\n', median(SCORE.Ea(:)));
    fprintf( '\t- Mean number of False Positives (over all voxels): %.2f\n', mean(SCORE.FP(:)));
    fprintf( '\t- Mean number of False Negatives (over all voxels): %.2f\n', mean(SCORE.FN(:)));
    fprintf( '\t- Mean Success Rate (over all voxels): %.2f\n', mean(SCORE.SUCCESS(:)));
    fprintf( '\t- Time ellapsed: %d seconds\n', time);
    fprintf( '\n');
    fprintf( '   ====================================================================== \n');

%% Save results

    if SAVE_RESULTS
     save(filename, 'FIBERS_rw', 'nFIBERS_rw','X', 'X_0','param', 'time', 'SCORE');

     save(globalfilename, 'SCORE');

     niiFIBERS.img=FIBERS;
     save_untouch_nii(niiFIBERS, niifname);
    end



