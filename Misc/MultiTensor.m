classdef MultiTensor

properties
	M		% number of fiber compartments
	f		% volume fractions
	R		% rotation matrices corresponding to each fiber orientation (theta,phi)
	lambda	% diffusivity profile for each fiber [mm^2/s]
end


methods

	% Constructor
	% -----------
	function obj = MultiTensor()
		obj.M			= 2;
		obj.f			= [.5 .5];
		obj.lambda		= [.3 .3 1.7 ; .3 .3 1.7]' * 1e-3;
		obj.R           = [];
		obj.R(:,:,1)	= obj.ROTATION(   0, pi/2);
		obj.R(:,:,2)	= obj.ROTATION(pi/2, pi/2);
	end


	% Compute the "TRUE ODF" corresponding to this model
	% --------------------------------------------------
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
% 		ODF = ODF / norm(ODF(:));
	end


	% estimate the FA for each fiber compartment
	% ------------------------------------------
	function fa = FA( obj )
		fa = zeros(1,obj.M);

		for i = 1:obj.M
			tmp = obj.lambda(:,i);
			l = sum( tmp(:) )/3;
			fa(i) = sqrt( 1.5 * sum((tmp(:)-l).^2) / sum(tmp(:).^2) );
		end
	end


	% Probe the signal at a given q-space coordinate
	% ----------------------------------------------
	function signal = E( obj, bCoord )
		if size(bCoord,1)==1, bCoord = bCoord'; end

		b = norm(bCoord);
		if b>0, bCoord = bCoord / b; end

		signal = 0;
		for i = 1:obj.M
			Di = obj.R(:,:,i) * diag(obj.lambda(:,i)) * obj.R(:,:,i)';
			signal = signal + obj.f(i) * exp(-b * bCoord' * Di * bCoord);
		end
	end


	% Add Rician noise to the signal
	% ------------------------------
	function En = addNoise( obj, E, sigma )
		if ( sigma<0 ), error('sigma must be >= 0'), end
		n = sigma * randn(1,2);
 		En = sqrt( (E+n(1))^2 + n(2)^2 );
	end


	% Probe the signal at several q-space positions
	% ---------------------------------------------
	function signal = acquireWithScheme( obj, grad_list, sigma )
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
			signal(d) = obj.E( XYZ(d,4) * XYZ(d,1:3)' );
			signal(d) = obj.addNoise( signal(d), sigma );
		end
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
