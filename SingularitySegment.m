classdef SingularitySegment
    % Defines a signularity segment of type (j), (j,+s), or (j,-s) 
    
    properties
        phi_i, phi_f
        th_i, th_f
    end
    
    methods
        function obj = j(obj, table, j)
            % IN : index j
            % OUT : segment (phi0,th0), (phif,thf)
            if (j < 0 || j > table.k) 
                error("j must be in [1, k]");
            end
            obj.phi_i = table.a(j);
            obj.phi_f = table.a(j);
            obj.th_i = 0;
            obj.th_f = pi;
        end

        function obj = j_s_neg(obj, table, j, s)
            if (j < 0 || j > table.k || s <= 0) 
                error("j must be in [1, k], or s > 0");
            end
            obj.phi_i = table.a(mod(j-2, table.k) + 1);  
            obj.phi_f = table.a(j);  
            if ( obj.phi_i > obj.phi_f)
                obj.phi_f = obj.phi_f + 2*pi;  % Periodicity of phi
            end
            obj.th_i = (obj.phi_f - obj.phi_i)./(2*s);
            obj.th_f = 0;

            eps = 1e-5;
            if (obj.th_i >= pi)
                obj.th_i = pi - eps;
                obj.phi_i = table.a(j) - 2*s*pi;
            end
        end

        function obj = j_s_pos(obj, table, j, s)
            if (j < 0 || j > table.k || s <= 0) 
                error("j must be in [1, k], or s > 0");
            end
            obj.phi_i = table.a(j);  
            obj.phi_f = table.a(mod(j, table.k) + 1);  
            if (obj.phi_i > obj.phi_f)
                obj.phi_f = obj.phi_f + 2*pi;  % Periodicity of phi
            end
            obj.th_i = 0;
            obj.th_f = (obj.phi_f - obj.phi_i)./(2*s);
            
            eps = 1e-6;
            if (obj.th_f >= pi)
                obj.th_f = pi - eps;
                obj.phi_f = table.a(j) + 2*s*pi - eps;
            end
        end

        function add_to_plot(obj, color)
            N = 10;
            x = mod(linspace(obj.phi_i, obj.phi_f, N), 2*pi+1e-4);
            y = linspace(obj.th_i, obj.th_f, N);
            if nargin < 2
                plot(x, y, 'LineWidth', 2);
            else
                plot(x, y, 'LineWidth', 2, 'Color', color);
            end
        end

        function [dom, im] = plot_image(obj, table, N, resolution, color)
            dom = linspace_vec([obj.phi_i, obj.th_i], [obj.phi_f, obj.th_f], resolution)';
                
            im = zeros(size(dom)); 
            for ii = 1:resolution
                orb = Orbit(table, dom(ii,1), dom(ii,2), N);
                im(ii,:) = orb.iter(end, :);
            end

            if nargin == 4
                plot(im(:,1), im(:,2), ".");
            else
                plot(im(:,1), im(:,2), ".", "Color", color);
            end
        end
        
        %{
        function [x, y] = intersection(obj, L2)
            dx1 = obj.phi_f - obj.phi_i;
            dy1 = obj.th_f - obj.th_i;
            dx2 = L2.phi_f - L2.phi_i;
            dy2 = L2.th_f - L2.th_i;
            
            A = [-dx1, dx2; -dy1, dy2];
            b = [L2.phi_i - obj.phi_i; L2.th_i - obj.th_i];
            if (rank(A) < 2)
                error("Cannot find intersection between segments!");
            end

            t = A \ b;
            x = obj.phi_i + t * dx1;
            y = obj.th_i + t * dy1;
        end
        %}
    end

    methods (Static)
        function [phi, th] = intersection_segments(table, j, s, t)
            % Lemma 15
            dj = table.b(j) - table.a(j);
            if ( s < 0 || t < 0 || s + t < dj/(2*pi))
                error("Intersection does not exist");
            end
            
            phi = table.a(j) + s*dj/(s+t);
            th = dj/(2*s + 2*t);
        end
    end
end

function arr = linspace_vec(vi, vf, N)
    arr = bsxfun(@plus,((vf(:)-vi(:))./(N-1))*[0:N-1],vi(:));
end