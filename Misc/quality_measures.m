function [SCORE] = quality_measures( FIBERS, FIBERS_gt, nFIBERS_gt, niiMASK )


    % [ SCORE ] = quality_measures( FIBERS, FIBERS_gt, nFIBERS_gt, niiMASK )
    % Computes quantitative measures to evaluate the global quality of the
    % reconstruction
    %
    %       INPUT:
    %       - FIBERS: reconstructed fiber configuration
    %       - FIBERS_gt: ground truth fiber configuration
    %       - nFIBERS_gt: 3D volume with the number of true fibers per voxel
    %       - MASK: indicates whether each voxel is active or not
    %
    %       OUTPUT:
    %         
    %       - SCORE:
    %           - Ea: (intra-voxel) Average Angular Error
    %           - SUCCESS: (intra-voxel) success(1)/fail(0) in identifying
    %           the fibres population
    %           - FP: (intra-voxel) n? of False Positives
    %           - FN: (intra-voxel) n? of False Negatives
    

% FIBERS and FIBERS_gt are two [(nX) x (nY) x (nZ) x (nK*3)] volumes. Each voxel contains up to nK peaks directions. 
% A peak direction is represented by a unit vector, scaled by the volume fraction of the corresponding peak.
nX=size(FIBERS,1);
nY=size(FIBERS,2);
nZ=size(FIBERS,3);

nK=size(FIBERS,4)/3; % nK: number of fibers per voxel // nK * 3 coordinates/fiber

nVOXELS=sum(niiMASK.img(:)); %number of (ACTIVE) voxels from which we want to recover the fiber configuration


% NUMBER OF FIBER POPULATIONS
% ===========================
M_true = nFIBERS_gt; %3D-vol with the number of true fibers in each voxel

MASK_fibers=zeros(size(M_true));
idf=find(M_true);
MASK_fibers(idf)=1; MASK_fibers=logical(MASK_fibers);

MASK=MASK_fibers & logical(niiMASK.img); 

SCORE.SUCCESS=zeros(size(M_true));
SCORE.FN=zeros(size(M_true));
SCORE.FP=zeros(size(M_true));
SCORE.Ea=zeros(size(M_true));


for ix=1:nX
    for iy=1:nY
        for iz=1:nZ
            
            if MASK(ix,iy,iz) == 0
                continue
            end
            
            % recovered Voxel configuration
            % =============================
            
            DIR=zeros(nK,3);
            for d=1:nK
                DIR(d,1:3)=FIBERS(ix,iy,iz,3*(d-1)+1:3*(d-1)+3);
                if ~any(DIR(d,:))
                    d=d-1;
                    break;
                end
            end
            DIR=DIR(1:d,:); % remove the "zero-directions"
            
            ESTIMATE = Voxel();
            ESTIMATE.M = size(DIR,1);
            for i=1:size(DIR,1)
                x=DIR(i,1); y=DIR(i,2); z=DIR(i,3);
                r= sqrt(x^2+y^2+z^2);
                theta= acos(z/r);
                phi=atan2(y,x);
                ESTIMATE.R(:,:,i) = ESTIMATE.ROTATION(phi,theta);
            end


            % true Voxel configuration
            % ==============================================================

            DIR_gt=zeros(nK,3);
            for d=1:nK
                DIR_gt(d,1:3)=FIBERS_gt(ix,iy,iz,3*(d-1)+1:3*(d-1)+3);
                if ~any(DIR_gt(d,:))
                    d=d-1;
                    break;
                end
            end
            DIR_gt=DIR_gt(1:d,:); % remove the "zero-directions"
            
            MODEL = Voxel();
            MODEL.M = size(DIR_gt,1);
            for i=1:size(DIR_gt,1)
                x=DIR_gt(i,1); y=DIR_gt(i,2); z=DIR_gt(i,3);
                r= sqrt(x^2+y^2+z^2);
                theta= acos(z/r);
                phi=atan2(y,x);
                MODEL.R(:,:,i) = MODEL.ROTATION(phi,theta);
            end 

            


            % ANGULAR RESOLUTION
            % ==================
            
            % precompute matrix with angular errors among all estimated and true fibers
            if (MODEL.M>0 && ESTIMATE.M>0)
                
                A = zeros( MODEL.M, ESTIMATE.M );
                for i = 1:MODEL.M
                    DIR_real = MODEL.R(:,:,i) * [0;0;1];
                    DIR_real = DIR_real / norm(DIR_real);
                    for j = 1:ESTIMATE.M
                        DIR_est = ESTIMATE.R(:,:,j) * [0;0;1];
                        DIR_est = DIR_est / norm(DIR_est);
                        err = acos( DIR_est' * DIR_real );
                        err = min( err, pi-err) / pi * 180;
                        A(i,j) = err;
                    end
                end
                

                % compute the "base" error
                Atmp = A;
                M = min(MODEL.M,ESTIMATE.M);
                err = zeros( 1, M );
                dT = 1:MODEL.M;
                dE = 1:ESTIMATE.M;
                for i = 1:M
                    [err(i),idx] = min( Atmp(:) );
                    [r c] = ind2sub( [MODEL.M ESTIMATE.M], idx );
                    Atmp(r,:) = NaN;
                    Atmp(:,c) = NaN;
                    dT( dT==r ) = [];
                    dE( dE==c ) = [];
                end

                % no penalisation for under/over-estimates 
                Ea = mean( err );

                SCORE.Ea(ix,iy,iz)=Ea;
                
                % compute FP and FN
                % -----------------
                TrueFiberAngleThr = 20; % acceptance cone [degree]

                SCORE.FP(ix,iy,iz) = 0; % FP - False Positives
                for j = 1:ESTIMATE.M
                   if nnz( A(:,j) < TrueFiberAngleThr ) == 0
                       SCORE.FP(ix,iy,iz) = SCORE.FP(ix,iy,iz) + 1;
                   end
                end

                SCORE.FN(ix,iy,iz) = 0; % FN - False Negatives
                for i = 1:MODEL.M
                   if nnz( A(i,:) < TrueFiberAngleThr ) == 0
                       SCORE.FN(ix,iy,iz) = SCORE.FN(ix,iy,iz) + 1;
                   end
                end
                
                SCORE.SUCCESS(ix,iy,iz)=0;
                if SCORE.FP(ix,iy,iz) == 0 && SCORE.FN(ix,iy,iz) == 0
                    SCORE.SUCCESS(ix,iy,iz)=1;
                end

                
            elseif (MODEL.M==0 && ESTIMATE.M>0)
                SCORE.SUCCESS(ix,iy,iz)=0;
                SCORE.FP(ix,iy,iz)=ESTIMATE.M;
                
            elseif (MODEL.M>0 && ESTIMATE.M==0)
                SCORE.SUCCESS(ix,iy,iz)=0;
                SCORE.FN(ix,iy,iz)=MODEL.M;
                
            end
            
            
            
        
        end
    end
end



end

