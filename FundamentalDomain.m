classdef FundamentalDomain
    %FUNDAMENTALDOMAIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vertices = zeros(3,2);  % Coordinates of vertices
        T = zeros(3,1);  % Connectivity
    end
    
    methods
        function obj = FundamentalDomain(table, j)
            L0 = SingularitySegment;
            L0 = L0.j(table, j);
            L1 = SingularitySegment;
            L1 = L1.j_s_pos(table, j, 1);

            obj.vertices(1,1) = L0.phi_i;
            obj.vertices(1,2) = L0.th_i;

            obj.vertices(2,1) = L1.phi_f;
            obj.vertices(2,2) = L1.th_f;
            
            obj.vertices(3,1) = L0.phi_i;
            obj.vertices(3,2) = L1.th_f;

            obj.T(1) = 2;
            obj.T(2) = 3;
            obj.T(3) = 1;
        end
        
        function add_to_plot(obj)
            hold on;
            pgon = polyshape(obj.vertices);
            plot(pgon);
        end

        function plot_sweep(obj, table, it, N)
            L = ExtendedSingularitySet(table);
            
            eps = 1e-6;
            top = linspace_vec(obj.vertices(2,:) + [eps, -eps], obj.vertices(3,:) + [eps, -eps], it)';
            Li = SingularitySegment;  % We will modify it but use its methods for plotting

            figure;
            figure(1); clf;
            L.new_plot; grid on; hold on;
            figure(2); clf;
            L.new_plot; grid on; hold on;
            
            colors = parula(it);
            
            resolution = 1000;
            svec = 1/2 * (top(:,1)-obj.vertices(1,1)) ./ (top(:,2)-obj.vertices(1,2));
            j = table.determine_arc(top(2,1));
            rj = table.r(j);
            rnext = table.r(mod(j,table.k) + 1);
            muj = sqrt(rj/rnext);
            for ii = 1:it
                Li.phi_i = obj.vertices(1,1);
                Li.th_i = obj.vertices(1,2);
                Li.phi_f = top(ii,1);
                Li.th_f = top(ii,2);

                figure(1); hold on;
                Li.add_to_plot(colors(ii,:)); 
                figure(2); hold on;
                [dom, im] = Li.plot_image(table, N, resolution, colors(ii, :));
                %figure(3); hold on;
                %subplot(2,1,1); hold on;
                %[m, id] = max(abs(im(:,2)-muj*dom(:,2)));
                %plot(svec(ii), m, ".", "MarkerSize", 14, "Color", colors(ii,:));  hold on;
                %subplot(2,1,2); hold on;
                %plot(svec(ii), table.determine_arc(im(id, 1)), ".", "MarkerSize", 14, "Color", colors(ii,:));
                ii = ii+1;
            end
            
            %figure(3);
            %subplot(2,1,1); grid on; ylim("Padded");
            %title("Maximum absolute change in theta with pre-image");
            %grid on; ylim("Padded"); xlim([0,1]); xlabel("corresponding s of domain");
            %subplot(2,1,2); grid on; ylim("Padded");
            %yticks([1,2,3,4])
            %title("Arc index of phi having the most change in theta");
        end

        function out = is_inside(obj, th, phi)
            % Uses barycentric coordinates
            x = th; y = phi;
            x1 = obj.vertices(1, 1); y1 = obj.vertices(1, 2);
            x2 = obj.vertices(2, 1); y2 = obj.vertices(2, 2);
            x3 = obj.vertices(3, 1); y3 = obj.vertices(3, 2);
        
            denom = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);
            alpha = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / denom;
            beta = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / denom;
            gamma = 1 - alpha - beta;
        
            out = (alpha >= 0) && (beta >= 0) && (gamma >= 0);
        end
    end
end


function arr = linspace_vec(vi, vf, N)
    arr = bsxfun(@plus,((vf(:)-vi(:))./(N-1))*[0:N-1],vi(:));
end