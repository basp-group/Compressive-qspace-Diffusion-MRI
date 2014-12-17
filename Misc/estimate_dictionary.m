clear variables, close all, fclose('all'); clc

%addpath(genpath('/Users/auria/Documents/MATLAB/spgl1-1.7'))
%addpath(genpath('/Users/auria/Documents/MATLAB/spams-matlab'))
%addpath(genpath('/Users/auria/Documents/MATLAB/FISTA_LASSO'))
addpath(genpath('Data')); %HARDI data (Challenge12)
addpath(genpath('NIFTI')); %nifti functions


addpath(genpath('Misc')); %miscellanious functions
addpath(genpath('Solvers')); %solvers
addpath(genpath('Operators')); %proximity operators

signal_path='Data/Measurements/';
data_path='Data/';
estimdic_path='Data/Estimated_Dicos/';
dtifit_path='Data/Estimated_Dicos/DTI_fit';

SNR			= 10;
nSAMPLES = 15;

fsignal=fullfile(signal_path,sprintf('ISBI12__signal_%d_%ddB.nii',nSAMPLES, SNR));
fout='Data/Estimated_Dicos/DTI_fit/dti';
fbvecs=sprintf('Data/bvecs_%d.txt', nSAMPLES);
fbvals=sprintf('Data/bvals_%d.txt', nSAMPLES);

fmask=fullfile(data_path,'ISBI12__mask.nii');

[status, result] = system( sprintf('. /usr/local/fsl/etc/fslconf/fsl.sh && export FSLOUTPUTTYPE=NIFTI && /usr/local/fsl/bin/dtifit --data=%s --out=%s --bvecs=%s --bvals=%s --mask=%s',fsignal,fout,fbvecs,fbvals,fmask));



niiFA = load_untouch_nii( fullfile(dtifit_path,'dti_FA.nii' ));
niiL1 = load_untouch_nii( fullfile(dtifit_path,'dti_L1.nii' ));
niiL2 = load_untouch_nii( fullfile(dtifit_path,'dti_L2.nii' ));
niiL3 = load_untouch_nii( fullfile(dtifit_path,'dti_L3.nii' ));
niiMASK = load_untouch_nii(fmask);
BVECS		= importdata( fullfile(fbvecs) );
BVALS		= importdata( fullfile(fbvals) ); 


IDX_WM = niiMASK.img==1 & niiFA.img>0.75;
V_WM = [ niiL1.img(IDX_WM) niiL2.img(IDX_WM) niiL3.img(IDX_WM) ];


KERNEL_WM = mean( V_WM, 1 );


if (1)
    figure(1), clf, hold on
    plot3( V_WM(:,1), V_WM(:,2), V_WM(:,3), 'x' )
    plot3( KERNEL_WM(1), KERNEL_WM(2), KERNEL_WM(3), 'ro' )
    axis equal
    xlabel('L1'), ylabel('L2'), zlabel('L3')
end


l=(KERNEL_WM(2)+KERNEL_WM(3))/2;
KERNEL_WM(3)=KERNEL_WM(1);
KERNEL_WM(2)=l;
KERNEL_WM(1)=l;


fprintf('KERNEL WM TO USE = [ %.2f; %.2f; %.2f ]\n', KERNEL_WM(3)*1e4, KERNEL_WM(2)*1e4, KERNEL_WM(1)*1e4)


%% create the DICTIONARY
%  =====================

DICTIONARY = [];
DICTIONARY.b	= BVALS(2:end);

nDIRECTIONS		= 200;
DIRECTIONS		= SphereGeneric( 'FromFile', 'gendir/gendir_200.txt', true );
ODF_SAMPLING	= SphereGeneric( 'golden', 724 );

XYZB=BVECS';XYZB=XYZB(2:end,:); %xyzb size (nSAMPLES,4)
XYZB(:,4) = BVALS(2:end);

fprintf( '-> Generating the overcomplete dictionary:\n' );

VOXEL   = Voxel();
VOXEL.M = 1;
VOXEL.f = 1;
VOXEL.lambda = KERNEL_WM';

for d = 1:nDIRECTIONS
	% rotate the template along the i-th direction
	VOXEL.R = VOXEL.ROTATION( DIRECTIONS.phi(d), DIRECTIONS.theta(d) );

	% simulate data acquisition of this configuration
	DICTIONARY.ATOMs(:,d) = VOXEL.acquireWithScheme( XYZB, 0, 'MultiTensor' );

	% generate the ODF for this configuration
	O = VOXEL.ODF_true( ODF_SAMPLING.x, ODF_SAMPLING.y, ODF_SAMPLING.z );
	DICTIONARY.ODFs(:,d)  = O(:);
end

% add the Isotropic kernel to the dictionary
d=d+1;
VOXEL.lambda = KERNEL_WM(3);
% rotate the template along the i-th direction
VOXEL.R = VOXEL.ROTATION( 0,0 );

% simulate data acquisition of this configuration
DICTIONARY.ATOMs(:,d) = VOXEL.acquireWithScheme( XYZB, 0, 'MultiTensor' );

% generate the ODF for this configuration
O = VOXEL.ODF_true( ODF_SAMPLING.x, ODF_SAMPLING.y, ODF_SAMPLING.z );
DICTIONARY.ODFs(:,d)  = O(:);
    
fprintf( '   [ OK ]\n' );


save(fullfile(estimdic_path,sprintf('Demo_DICTIONARY_%d_%d_%d',nDIRECTIONS, nSAMPLES, SNR)),'DICTIONARY');

