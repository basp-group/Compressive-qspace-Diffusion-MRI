classdef SphereGeneric < handle
	properties
        type, L
		x, y, z, tri
		theta, phi
	end
	
	methods
		%% constructor
		function obj = SphereGeneric( gridType, n, isHalf )
			if isequal(gridType,'golden')
				if mod(n,2) || n<10
					error('"n" parameter must be EVEN and > 9')
				end
				inc = pi * (3-sqrt(5));
				off = 2/n;
				k = 0:n-1;
				y = k * off - 1 + (off/2);
				r = sqrt(1 - (y.^2));
				phi = k * inc;
				p = [cos(phi).*r; y ;sin(phi).*r]';

				for d = 1:size(p,1)/2
					p(d,:) = p(d,:) / norm(p(d,:));
                    obj.x(d,1) = p(d,1);
                    obj.y(d,1) = p(d,2);
                    obj.z(d,1) = p(d,3);
                    obj.x(d+size(p,1)/2,1) = -p(d,1);
                    obj.y(d+size(p,1)/2,1) = -p(d,2);
                    obj.z(d+size(p,1)/2,1) = -p(d,3);
					
					r = sqrt( p(d,1)^2 + p(d,2)^2 + p(d,3)^2 );
					obj.theta(d,1) = acos(p(d,3)/r);
					obj.theta(d+size(p,1)/2,1) = obj.theta(d,1);
					obj.phi(d,1)   = atan2(p(d,2),p(d,1));
					obj.phi(d+size(p,1)/2,1)   = obj.phi(d,1);
				end
				
				obj.tri = convhulln( [obj.x obj.y obj.z] );
                
                % anna ==========================
                if nargin<3, isHalf = false; end
				if isHalf
					obj.x = [ obj.x ; -obj.x ];
					obj.y = [ obj.y ; -obj.y ];
					obj.z = [ obj.z ; -obj.z ];
                    [obj.theta, obj.phi] = ssht_c2s( obj.x, obj.y, obj.z );
                end
                % =============================

			elseif  isequal(gridType,'FromFile')
				% "n" is used here as the filename
				SCRIPTPATH = fileparts( mfilename('fullpath') );
				if exist( n )
					G = importdata( n );
				else
					G = importdata( fullfile(SCRIPTPATH,n) );
				end
				obj.x = G(:,1);
				obj.y = G(:,2);
				obj.z = G(:,3);
	
				if nargin<3, isHalf = false; end
				if isHalf
					obj.x = [ obj.x ; -obj.x ];
					obj.y = [ obj.y ; -obj.y ];
					obj.z = [ obj.z ; -obj.z ];
				end
				
				[obj.theta, obj.phi] = ssht_c2s( obj.x, obj.y, obj.z );
				obj.tri = convhulln( [obj.x obj.y obj.z] );
				n = size(obj.x,1);
            
            elseif  isequal(gridType,'FromMatrix')
                % "n" is used here as the matrix name
                obj.x = n(:,1);
				obj.y = n(:,2);
				obj.z = n(:,3);
                
                if nargin<3, isHalf = false; end
				if isHalf
					obj.x = [ obj.x ; -obj.x ];
					obj.y = [ obj.y ; -obj.y ];
					obj.z = [ obj.z ; -obj.z ];
				end
				
				[obj.theta, obj.phi] = ssht_c2s( obj.x, obj.y, obj.z );
				obj.tri = convhulln( [obj.x obj.y obj.z] );
				n = size(obj.x,1);
            
			else
				error('"gridType" parameter not valid')
			end
			
			obj.type = gridType;
			obj.L    = n;
		end


		%% Plot this function
		%     passing MODEL also plot the fibers' information
		function plot( obj, values, plotTYPE, MODEL )
			if nargin<3, plotTYPE=1; end
			if nargin<4, MODEL=[]; end
            
            if ~isequal(size(obj.x),size(values))
                error('SphereGeneric.plot(): values''s size doesn''t match sphere''s size')
            end
			
			if ( plotTYPE==1 )
				% SPHERE visualization
				trisurf(obj.tri, obj.x,obj.y,obj.z, values)
				maxLEN = 1;
			elseif ( plotTYPE==2 )
                % ODF visualization
				l = repmat(values,1,3) .* [obj.x,obj.y,obj.z];
				maxLEN = max(sqrt( sum(l.^2,2) ));
				trisurf(obj.tri, l(:,1),l(:,2),l(:,3), values)
			elseif ( plotTYPE==3 )
				% ODF visualization with min-max normalization
				MIN = min(values);
				MAX = max(values);
 				l = repmat((values-MIN) / (MAX-MIN),1,3) .* [obj.x,obj.y,obj.z];
				maxLEN = max(sqrt( sum(l.^2,2) ));
				trisurf(obj.tri, l(:,1),l(:,2),l(:,3), values)
			elseif ( plotTYPE==0 )
				maxLEN = 1;
			else
				return
			end

			
			% plot model's fiber directions
			if ~isempty(MODEL)
				hold on
				h = zeros(1,MODEL.M);
				LEGEND = cell(1,MODEL.M);
				for i = 1:MODEL.M
					xyz = MODEL.R(:,:,i) * [0;0;1] *  maxLEN * 1.2;
					h(i) = line( [-xyz(1) xyz(1)], [-xyz(2) xyz(2)], [-xyz(3) xyz(3)], 'Color',[0 0 0], 'LineWidth',2 );
					LEGEND{i} = sprintf('f=%.2f',MODEL.f(i));
				end
				legend( h, LEGEND ); legend('off');
			end
			
			view(3), rotate3d on, axis equal, shading interp; lighting phong
			xlabel('x'), ylabel('y'), zlabel('z')
			colormap(jet(256));
			
			drawnow
		end

	end
end
