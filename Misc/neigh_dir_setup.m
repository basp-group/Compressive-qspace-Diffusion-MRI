function [ NEIGH, DISTmin ] = neigh_dir_setup( DIRECTIONS, DEGREE_max )

    % [ NEIGH, DISTmin ] = neigh_dir_setup( DIRECTIONS, DEGREE_max )
    % creates a cell array containing, in each cell -indexed I-, the vector of indices corresponding to the
    % neighbor-directions (within an angular range of DEGREE_min) of direction I
    %
    %   INPUTS:
    %
    %       - DIRECTIONS: SphereGeneric object storing a set of directions
    %       on the sphere
    %
    %       - DEGREE_max: angular threshold to define the neighborhood (in degrees) 
    %
    %   OUTPUT:
    %
    %       - NEIGH: cell array indicating the nearest
    %       neighbor directions for each direction in object DIRECTIONS. 
    %       


ANGLE_thr	= DEGREE_max/180*pi;
nDIRECTIONS = numel(DIRECTIONS.theta);
NEIGH		= cell(nDIRECTIONS/2,1);
DISTmin=pi;

for d1 = 1:nDIRECTIONS/2
	DIST = zeros(nDIRECTIONS,1);
	for d2 = 1:nDIRECTIONS
		a = [DIRECTIONS.x(d1) DIRECTIONS.y(d1) DIRECTIONS.z(d1)];
		b = [DIRECTIONS.x(d2) DIRECTIONS.y(d2) DIRECTIONS.z(d2)];
		DIST(d2) = atan2(norm(cross(a,b)),dot(a,b));
        
        if (DIST(d2)>pi/2) DIST(d2)=pi-DIST(d2); end %the direction of the vector is not relevant
        if (DIST(d2)<DISTmin && DIST(d2)>0), DISTmin=DIST(d2); end
	end
	NEIGH{d1} = find( DIST>0 & DIST <= ANGLE_thr );
    NEIGH{d1} = mod(NEIGH{d1},nDIRECTIONS/2); 
    NEIGH{d1}(NEIGH{d1}==0)=nDIRECTIONS/2; 
    
end

