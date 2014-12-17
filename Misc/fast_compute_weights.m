function [ W ] = fast_compute_weights( X, NEIGH_DIR, NEIGH_VOX, MASK ) 

% M-function wrapper for parsing the inputs



    if nargin < 4
        if nargin <3
            error('Not enough input arguments');
        end
        MASK=ones(size(X,2)); MASK=logical(MASK); %tau is unspecified, set default value (mask open for all voxels)
    end
    


if ndims(X)~=2 || ~isreal(X) || issparse(X) || ~isa(X,'double')
    error('X must be a real 2D full double array');
elseif ~isa(NEIGH_DIR,'cell') || numel(NEIGH_DIR)~=size(X,1)
    error('NEIGH_DIR must be a cell array with the same dimension as the number of rows of X');
elseif ~isa(NEIGH_VOX,'cell') || numel(NEIGH_VOX)~=size(X,2)
    error('NEIGH_VOX must be a cell array with the same dimension as the number of columns of X');
elseif ~islogical(MASK) || numel(MASK)~=size(X,2)
    error('MASK must be a logical array with the same number of elements as the number of columns of X');

end

MASK=MASK(:);

mex('compute_weights_mex.c');
W=compute_weights_mex( X, NEIGH_DIR, NEIGH_VOX, MASK );
        
end

