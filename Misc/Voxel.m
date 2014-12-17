classdef Voxel

properties
	M				% number of fiber compartments
	f				% volume fractions
	R				% rotation matrices corresponding to each fiber orientation (theta,phi)

	% used if 'MultiTensor'
	lambda			% diffusivity profile for each fiber [mm^2/s]
	
	% used if 'Soderman'
	r, l			% radius and length of the cylinder(s)
	D				% diffusivity of the molecules
	Delta, delta	% acquisition parameters
end


methods

	% Constructor
	% -----------
	function obj = Voxel()
		obj.M			= 2;
		obj.f			= [.5 .5];
		obj.R           = [];
		obj.R(:,:,1)	= obj.ROTATION(   0, pi/2);
		obj.R(:,:,2)	= obj.ROTATION(pi/2, pi/2);
		
		% parameters for Soderman's model (as in Ozarslan 2006)
		obj.r			= 5e-3;
		obj.l			= 5;
		obj.D			= 2.02e-3;
		obj.Delta		= 20.8e-3;
		obj.delta		= 2.4e-3;
		
		% tensor corresponding to the previous parameters (estimated from data)
		obj.lambda		= [3.194e-04 3.194e-04 2.093e-03; 3.194e-04 3.194e-04 2.093e-03]';
	end



	% Probe the signal at several q-space positions
	% ---------------------------------------------
	function signal = acquireWithScheme( obj, grad_list, sigma, model )
		if nargin < 4
			model = 'MultiTensor';
		end

		if isstr(grad_list)
			if ~exist(grad_list,'file')
				error('unable to locate file');
			end
			XYZ = importdata( grad_list, ' ' );
		else
			if size(grad_list,2) ~= 4
				error('the gradient list must be (nx4)');
			end
			XYZ = grad_list;
		end

		nDIR   = size( XYZ, 1 );
		signal = zeros( nDIR, 1 );

		for d = 1:nDIR
			signal(d) = obj.E( XYZ(d,4) * XYZ(d,1:3)', model );
			signal(d) = obj.addNoise( signal(d), sigma );
		end
	end



	% Probe the signal at a given q-space coordinate
	% ----------------------------------------------
	function signal = E( obj, bCoord, model )
		if nargin < 3
			model = 'MultiTensor';
		end
		if size(bCoord,1)==1, bCoord = bCoord'; end
		
		b = norm(bCoord);
		if b<1,	signal = 1;	return;  end
		if b>0, bCoord = bCoord / b; end

		signal = 0;
		switch ( model )
			
			case 'MultiTensor'		
				for i = 1:obj.M
					Di = obj.R(:,:,i) * diag(obj.lambda(:,i)) * obj.R(:,:,i)';
					signal = signal + obj.f(i) * exp(-b * bCoord' * Di * bCoord);
				end

			case 'Soderman'
				q = sqrt( b / ( 4.0*pi^2 * (obj.Delta-obj.delta/3.0) ) );
				for i = 1:obj.M
					dir = obj.R(:,:,i) * [0;0;1];
					dir = dir ./ norm(dir);
					theta = acos( abs( dir' * bCoord ) );
					if theta>pi/2
						theta = mod( theta, pi/2 );
					end
					if abs(theta) < 1e-9, theta = 1e-9; end
					if abs(theta-pi/2) < 1e-9, theta = pi/2 - 1e-9; end

					% derivatives of Bessel functions
					x = 2*pi * q *obj.r * sin(theta);
					Jp = zeros( 10, 1 );
					for m = 0:9
						Jp(m+1) = ( 0.5*( besselj(m-1,x) - besselj(m+1,x) ) ) .^ 2;
					end

					signal = signal + obj.f(i) * Voxel_soderman( q, theta, obj.r, obj.l, obj.Delta, obj.D, Jp );
				end

			otherwise
				error( 'Voxel::E(): impossible to simulate signal for this model' )
		end
	end



	% Add Rician noise to the signal
	% ------------------------------
	function En = addNoise( obj, E, sigma )
		if ( sigma<0 ), error('Voxel::addNoise(): sigma must be >= 0'), end
		n = sigma * randn(1,2);
 		En = sqrt( (E+n(1))^2 + n(2)^2 );
	end

	

	% Compute the "TRUE ODF" corresponding to 'MultiTensor' model
	% -----------------------------------------------------------
	function ODF = ODF_true( obj, x, y, z )
		nDirections = size(x,1);
		ODF = zeros( nDirections, 1 );
		for idx = 1:nDirections
			r = [x(idx); y(idx); z(idx)];
			ODF(idx) = 0;
			for i = 1:obj.M
				Di = obj.R(:,:,i) * diag(obj.lambda(:,i)) * obj.R(:,:,i)';
				ODF(idx) = ODF(idx) + obj.f(i) / (4*pi*sqrt(det(Di))) * 1/sqrt( r' * inv(Di) * r ).^3;
			end
		end
 		ODF = ODF / sum(ODF(:));
	end



	% Compute a rotation matrix corresponding to the orientation (azimuth,zenith)
	%     azimuth (phi):	angle in the x-y plane
	%     zenith  (theta):	angle from z axis
	% ---------------------------------------------------------------------------
	function M = ROTATION( obj, azimuth, zenith )
		azimuth = mod(azimuth,2*pi);
		zenith  = mod(zenith,pi);
		M = [ cos(azimuth) -sin(azimuth) 0 ; sin(azimuth) cos(azimuth) 0 ; 0 0 1 ] * ...
		    [ cos(zenith) 0 sin(zenith) ; 0 1 0 ; -sin(zenith) 0 cos(zenith) ];
	end
end

end
