function [ act_NEIGH_VOX ] = map_neigh_vox( NEIGH_VOX, allV_id, actV_id )
%
%   map_neigh_vox "translates" the neighboring map NEIGH_VOX (whole 3D-cube
%   as a system of reference) to act_NEIGH_VOX, in which the indexes
%   correspond to the active VOXELS reference system
%
%   INPUTS:
%
%       -NEIGH_VOX: cell array (of dimension nVOXELS) indicating the nearest
%       neighbor voxels for each voxel. Indices references to the nVOXELS
%       3D-vol
%
%       -allV_id: vector (of dimension actVOXELS), for each active
%       voxel it contains its index in the global nVOXELS 3D-vol.
%
%       -actV_id: 3d-vol (of dimension nVOXELS), for each voxel, it
%       contains its index in the actVOXELS system of reference (in case it
%       is an active voxel) or zero (in case it is a non active voxel).
%
%   OUTPUT:
%
%       - act_NEIGH_VOX: cell array (of dimension actVOXELS) indicating the nearest
%       neighbor voxels for each voxel. Indices references to the actVOXELS
%       


actVOXELS=numel(allV_id);

act_NEIGH_VOX=cell(actVOXELS, 1);

for i=1:actVOXELS
    
    
    neigh=NEIGH_VOX{allV_id(i)};
    
    for j=1:numel(neigh)
        
        if actV_id(neigh(j)), act_NEIGH_VOX{i}=[act_NEIGH_VOX{i};actV_id(neigh(j))]; end
    end
    
    if isempty(find(act_NEIGH_VOX{i}==i,1)), act_NEIGH_VOX{i}=[act_NEIGH_VOX{i};i]; end
    
end
        
    
    



end

