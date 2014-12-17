function [ NEIGH ] = neigh_vox_setup( VOXELS, r )

    % [ NEIGH ] = neigh_vox_setup( VOXELS, r )
    % creates a cell array containing, in each cell -indexed V-, the vector of indices corresponding to the
    % neighboring voxels (within a radius r-voxels) of voxel V
    %
    %   INPUTS:
    %
    %       - VOXELS: 3D-array with the corresponding dimensions
    %
    %       - r: to define the neighborhood  
    %
    %   OUTPUT:
    %
    %       - NEIGH: cell array indicating the nearest
    %       neighbor voxels for each voxel. 
    %       



nVOXELS = numel(VOXELS(:));
dimX=size(VOXELS,1);
dimY=size(VOXELS,2);
dimZ=size(VOXELS,3);
NEIGH		= cell(nVOXELS,1);


for n=1:nVOXELS
    
    [i,j,k]=ind2sub(size(VOXELS),n);
    
    for dx=-r:r
        for dy=-r:r
            for dz=-r:r
                if (i+dx)>0 && (j+dy)>0 && (k+dz)>0 &&(i+dx)<=dimX && (j+dy)<=dimY && (k+dz)<=dimZ % && ~isequal([dx,dy,dz],[0,0,0])
                    NEIGH{n}=[NEIGH{n};sub2ind(size(VOXELS),i+dx,j+dy,k+dz)];
                end
            end
        end
    end
    
    
    
end

end

