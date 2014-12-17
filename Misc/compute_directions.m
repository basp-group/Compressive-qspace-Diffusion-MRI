function [ FIBERS, nFIBERS ] = compute_directions( X, niiMASK,NEIGH_DIR, DIRECTIONS,niiK )

        %
        % [ FIBERS, nFIBERS ] = compute_directions( X, niiMASK,NEIGH_DIR, DIRECTIONS,nDIRECTIONS, niiK )
        %  computes the main direction of the fiber(s) for every active voxel (mask open)
        %
        %       INPUT
        %       - X: matrix of dimensions (nATOMS, nVOXELS) containing the FOD for every voxel
        %       - niiMASK: mask
        %       - NEIGH_DIR: cell array of dimensions (nDIRECTIONS,1), for each
        %       direction, it stores the vector of indices corresponding to its
        %       nearest neighbors
        %       - DIRECTIONS: SphereGeneric object storing the directions of the
        %       atoms of our dictionary
        %       - niiK: number of expected fibers per voxel (nii format)
        %
        %       OUTPUT
        %       - FIBERS: 4-D array 
        %       - nFIBERS: 3-D volume indicating the number of recovered fibers per voxel


nX=niiMASK.hdr.dime.dim(2);
nY=niiMASK.hdr.dime.dim(3);
nZ=niiMASK.hdr.dime.dim(4);

FIBERS=zeros(niiMASK.hdr.dime.dim(2),niiMASK.hdr.dime.dim(3),niiMASK.hdr.dime.dim(4), niiK*3); %all fiber directions per voxel 
nFIBERS=zeros(niiMASK.hdr.dime.dim(2),niiMASK.hdr.dime.dim(3),niiMASK.hdr.dime.dim(4)); %number of fibers per voxel

res_ang=20;

for ix=1:nX
    for iy=1:nY
        for iz=1:nZ
            
            if niiMASK.img(ix,iy,iz) == 0
                 continue
            end

            v=sub2ind(size(rand(nX,nY,nZ)),ix,iy,iz);
            x=X(:,v); % x corresponds to the recovered FOD in voxel v
            
            [ DIR, MAX, idxMAX ] = maxima_detect( [x;x], DIRECTIONS, NEIGH_DIR, 1e-1 );
 
            % Save directions
            for d=1:size(DIR,1)

                if ~any(DIR(d,:))
                    d=d-1;
                    break;
                end
                FIBERS(ix,iy,iz,3*(d-1)+1:3*(d-1)+3)=DIR(d,1:3)*MAX(d);

            end

            if isempty(DIR), d=0; end
            if isempty(d), d=0; end

            nFIBERS(ix,iy,iz)=d;

            while (d<niiK)

                d=d+1;
                FIBERS(ix,iy,iz,3*(d-1)+1:3*(d-1)+3)=[0 0 0];

            end



            
        end
    end
end

end

