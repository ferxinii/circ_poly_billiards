classdef BilliardTable < handle
    
    properties
        k
        r = []
        a = []
        b = []
        O = []
    end
    
    methods
        function set_moss_egg(obj, r)
            l_r1 = r;  %l_ for "local"
            l_r2 = 2*r;
            l_r3 = (2-sqrt(2))*r;
            l_r4 = 2*r;
            d1 = pi;
            d2 = pi / 4;
            d3 = pi / 2;
            d4 = pi / 4;
            set_table_4(obj, l_r1, l_r2, l_r3, l_r4, d1, d2, d3, d4);
        end

        function s(obj, r, R, a)
            l_r1 = r;
            l_r3 = r;
            l_r2 = R;
            l_r4 = R;
            d1 = a;
            d3 = a;
            d2 = pi - a;
            d4 = pi - a;
            set_table_4(obj, l_r1, l_r2, l_r3, l_r4, d1, d2, d3, d4);
        end

        function set_table_4(obj, r1, r2, r3, r4, d1, d2, d3, d4)
            obj.k = 4;

            obj.r = zeros(4,1);
            obj.r(1) = r1;    obj.r(2) = r2;
            obj.r(3) = r3;    obj.r(4) = r4;
            
            obj.a = zeros(4,1);
            obj.a(1) = 0;               obj.a(2)= d1;
            obj.a(3) = obj.a(2) + d2;   obj.a(4) = obj.a(3) + d3;
            
            obj.b = zeros(4,1);
            obj.b(1) = d1;              obj.b(2) = obj.b(1) + d2;
            obj.b(3) = obj.b(2) + d3;   obj.b(4) = obj.b(3) + d4;
            
            obj.O = zeros(4,2);
            obj.O(1,1) = 0;              
            obj.O(1,2) = 0;
            obj.O(2,1) = obj.O(1,1) + (r1-r2)*cos(obj.b(1));
            obj.O(2,2) = obj.O(1,2) + (r1-r2)*sin(obj.b(1));
            obj.O(3,1) = obj.O(2,1) + (r2-r3)*cos(obj.b(2));
            obj.O(3,2) = obj.O(2,2) + (r2-r3)*sin(obj.b(2));
            obj.O(4,1) = obj.O(3,1) + (r3-r4)*cos(obj.b(3));
            obj.O(4,2) = obj.O(3,2) + (r3-r4)*sin(obj.b(3));
        end

        function set_table_6(obj, A, B, C, r)
            % A, B, C are the points that form a triangle
            al = norm(C-B);
            bl = norm(C-A);
            cl = norm(A-B);
            if (r <= max(max(0, al-cl), al-bl))
                error("r <= max(0, a-c, a-b)");
            end

            alpha = acos((bl^2 + cl^2 - al^2)/(2*bl*cl));
            beta = asin(bl/al*sin(alpha)); %acos((al^2 + cl^2 - bl^2)/(2*al*cl));
            gamma = asin(cl/al*sin(alpha)); %acos((al^2 + bl^2 - cl^2)/(2*al*bl));
            
            d1 = beta;      d4 = beta;
            d2 = gamma;     d5 = gamma;
            d3 = alpha;     d6 = alpha;

            obj.k = 6;
            obj.r = [r+cl; r+cl-al; r+cl-al+bl; r+bl-al; r+bl; r];

            obj.a = zeros(6,1);
            obj.a(1) = 0;               obj.a(2)= d1;
            obj.a(3) = obj.a(2) + d2;   obj.a(4) = obj.a(3) + d3;
            obj.a(5) = obj.a(4) + d4;   obj.a(6) = obj.a(5) + d5;
            
            obj.b = zeros(6,1);
            obj.b(1) = d1;              obj.b(2) = obj.b(1) + d2;
            obj.b(3) = obj.b(2) + d3;   obj.b(4) = obj.b(3) + d4;
            obj.b(5) = obj.b(4) + d5;   obj.b(6) = obj.b(5) + d6;

            obj.O = zeros(6,2);
            obj.O(1,:) = B;     obj.O(4,:) = B;
            obj.O(2,:) = C;     obj.O(5,:) = C;
            obj.O(3,:) = A;     obj.O(6,:) = A;
        end

        function [x, y] = polar_parametrization(obj, s_vec)
            % IN: vector of values, domain [0, 2pi]
            % OUT: vectors x(s), y(s)
            x = zeros(size(s_vec));
            y = zeros(size(s_vec));
            for ii = 1:length(s_vec)
                s = mod(s_vec(ii), 2*pi);
                j = obj.determine_arc(s);

                x(ii) = obj.O(j,1) + obj.r(j) * cos(s);
                y(ii) = obj.O(j,2) + obj.r(j) * sin(s);
            end
        end

        function [x, y] = polar_parametrization_deriv(obj, s_vec)
            % IN: vector of values, domain [0, 2pi]
            % OUT: vectors x'(s), y'(s)
            x = zeros(size(s_vec));
            y = zeros(size(s_vec));
            for ii = 1:length(s_vec)
                s = mod(s_vec(ii), 2*pi);
                j = obj.determine_arc(s);

                x(ii) = -obj.r(j) * sin(s);
                y(ii) = obj.r(j) * cos(s);
            end
        end
        
        function new_plot_detailed(obj)
            % Plots table in a new figure, detailing each arc
            figure();  hold on; fontsize(12,"points");
            grid on; axis equal; xlabel("x"); ylabel("y");
            xlim("padded"); ylim("padded");
            N = 1000;
            
            for ii = 1:obj.k
                % Arc, with appropriate number of points
                Ni = round(N*(obj.b(ii) - obj.a(ii))/(2*pi));
                s = linspace(obj.a(ii), obj.b(ii), Ni);
                [x, y] = polar_parametrization(obj, s);
                plot(x, y, 'LineWidth', 2, 'DisplayName', string(ii));
                
                % Singularities
                [x, y] = polar_parametrization(obj, obj.a(ii));
                plot(x, y, 'o', 'LineWidth', 2, 'Color', 'black', ...
                    'HandleVisibility','off');
            end
        end

        function add_to_plot(obj)
            N = 1000;
            s = linspace(0, 2*pi, N);
            [x, y] = polar_parametrization(obj, s);
            plot(x, y, 'LineWidth', 2, 'Color', 'black');
        end

        function add_to_plot_extended(obj)
            hold on;
            N = 100;
            s = linspace(0, 2*pi, N);
            for ii = 1:obj.k
                x = obj.O(ii,1) + obj.r(ii)*cos(s);
                y = obj.O(ii,2) + obj.r(ii)*sin(s);
                plot(x, y, 'LineWidth', 1);
            end
        end

        function j = determine_arc(obj, phi)
            if (phi >= obj.b(obj.k))
                error("phi >= bk, TODO?")
            end

            j = 0;
            for ii = 1:obj.k
                if (obj.is_in_arc(ii, phi))
                    j = ii;
                    break
                end
            end
            if (j == 0)
                error("cannot determine current arc. TODO?");
            end
        end

        function out = is_in_arc(obj, j, phi)
            eps = 1e-12;
            phi = mod(phi, 2*pi);
            cond1 = phi > obj.a(j) || abs(phi - obj.a(j)) < eps;
            cond2 = phi < obj.b(j);
            if (cond1  && cond2)
                out = 1;
            else 
                out = 0;
            end
        end

        function out = chi_min(obj, j)
            dj = obj.b(j) - obj.a(j);
            dnext = obj.b(mod(j, obj.k) + 1) - obj.a(mod(j, obj.k) + 1);
            rj = obj.r(j);
            rnext = obj.r(mod(j, obj.k) + 1);
            muj = sqrt(rj / rnext);

            out = 1 + ceil(muj * dj / (2 * dnext));
        end
    end
end

