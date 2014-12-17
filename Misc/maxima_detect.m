function [ DIR, MAX, idxMAX ] = maxima_detect( coeff, DIRECTIONS, NEIGH, THR )

if ( nargin<4 ), THR = 0.2; end

DIR    = [];
MAX    = [];
idxMAX = [];

f=coeff(1:numel(DIRECTIONS.theta));

nATOMS=numel(coeff)/2;

nDIRECTIONS = numel(DIRECTIONS.theta)/2;


for d1 = 1:nDIRECTIONS
	if all( f(d1) > f( NEIGH{d1} ) )
		idxMAX(end+1)  = d1;
		MAX(end+1)     = f(d1);
		DIR(end+1,1:3) = [ DIRECTIONS.x(d1) DIRECTIONS.y(d1) DIRECTIONS.z(d1) ];
	end
end

% Keep at most 5 directions
if numel(MAX)>5
    [sortVal,sortId] = sort(MAX(:),'descend');  %# Sort the values in descending order
    keep = sortId(1:5);  %# linear index into MAX of the 5 largest values

    MAX=MAX(keep);
    idxMAX=idxMAX(keep);
    DIR=DIR(keep,:);
end

% Keep only those which max val are at least THR*max(MAX);
keep   = ( MAX >= THR*max(MAX) );

DIR    = DIR( keep, : );
MAX    = MAX( keep );
idxMAX = idxMAX( keep );


