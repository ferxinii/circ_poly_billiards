classdef Orbit
    
    properties
        N  % Number of iterates
        iter = []  % iterates in phase space
    end
    
    methods
        function obj = Orbit(table, phi0, th0, N)
            inv = 0;
            if (N < 0)
                inv = 1;
                th0 = pi - th0;
                N = -N;
            end
            if (th0 > pi)
                error("theta > pi!");
            end
            
            obj.N = N + 1;
            obj.iter = zeros(N + 1, 2);

            phi = phi0;  th = th0;
            for ii = 1:obj.N
                obj.iter(ii, 1) = phi;
                obj.iter(ii, 2) = th;
                if inv == 1
                    obj.iter(ii, 2) = pi - th;
                end
                [phi, th] = obj.map(table, phi, th);
            end
        end

        function plot_cartesian(obj, table)
            hold on;
            for ii = 1:(obj.N - 1)
                [x0, y0] = table.polar_parametrization(obj.iter(ii, 1));
                [xf, yf] = table.polar_parametrization(obj.iter(ii+1, 1));
                plot([x0, xf], [y0, yf]);
            end
        end

        function plot_phasespace(obj, color)
            hold on
            if nargin == 2
                plot(obj.iter(:,1), obj.iter(:,2), '.', 'MarkerSize', 12, "Color", color);
            else
                plot(obj.iter(:,1), obj.iter(:,2), '.', 'MarkerSize', 12);
            end
        end

        function out = is_generic(obj, table)
            eps = 1e-12;
            out = 1;
            L = ExtendedSingularitySet(table);
            for ii = 1:obj.N
                phi = obj.iter(ii, 1);
                th = obj.iter(ii, 2);
                for jj = 1:length(L.L)
                    Lj = L.L(jj);
                    % Dimension z equal to 0, to use matlab's functions
                    ac = [Lj.phi_f - Lj.phi_i, Lj.th_f - Lj.th_i, 0];
                    ab = [phi - Lj.phi_i, th - Lj.th_i, 0];
                    
                    % Check colinearity and projection
                    aux1 = cross(ab, ac);
                    aux2 = dot(ab, ac);
                    if (all(abs(aux1) < eps) && aux2 > 0 && aux2 < norm(ac))
                        out = 0;
                        return;
                    end
                end
            end
        end

        function out = is_sliding(obj, table)
            out = 1;
            for ii = 1:obj.N-1
                current_arc = table.determine_arc(obj.iter(ii, 1));
                next_arc = table.determine_arc(obj.iter(ii+1, 1));
                diff = next_arc - current_arc;
                if (current_arc == table.k && next_arc == 1)
                    diff = 1;
                end
                if (diff ~= 0 && diff ~= 1)
                    out = 0;
                    return;
                end
            end
        end

        function out = extract_return_map(obj, table)
            D1 = FundamentalDomain(table, 1);
            D1.is_inside(obj.iter(1,1), obj.iter(1,2))
            if (obj.is_generic(table) == 0 || obj.is_sliding(table) == 0 || D1.is_inside(obj.iter(1,1), obj.iter(1,2)) == 0)
                error("Can only extract return map from fundamental generic sliding orbits")
            end
            
            arr = [];
            for ii = 1:obj.N
                if (D1.is_inside(obj.iter(ii,1), obj.iter(ii,2)) == 1)
                    arr = [arr; obj.iter(ii,1), obj.iter(ii,2)];
                end
            end

            out = Orbit(table, 0,0,1);
            out.N = length(arr);
            out.iter = arr;
        end
    end

    methods (Static)
        function [phi, th] = map(table, phi0, th0)
            if (th0 < 0 || th0 > pi)
                error("th0 must be in [0, pi]");
            end
            if (th0 == 0 || th0 == pi)
                phi = phi0;
                th = th0;
                return;
            end
            phi = Orbit.map_phi(table, phi0, th0);
            [xp, yp] = table.polar_parametrization_deriv(phi);
            ds = [xp, yp];

            j = table.determine_arc(phi0);
            vx = - table.r(j) * sin(phi0+th0);
            vy = table.r(j) * cos(phi0+th0);
            v = [vx, vy];

            arg = dot(ds, v)/(norm(ds)*norm(v));
            if (arg > 1 || arg < -1)
                error("Cannot compute th image");
            end

            th = acos(arg);
            while (th < 0)
                th = th + 2*pi
            end
            while (th > 2*pi)
                th = th - 2*pi
            end
            if (th > pi)
                th = 2*pi - th
            end
        end

        function out = map_phi(table, phi0, th0)
            if (phi0 >= 2*pi)
                error("phi0 must be in [0, 2pi)");
            end
            eps = 1e-12;
            j = table.determine_arc(phi0);
            % Check case colision is in node
            if (abs(phi0 + 2*th0 - table.b(j)) < eps)
                %disp("here");
                out = mod(phi0 + 2*th0, 2*pi);
                return;
            end

            [k_vec, phi1, phi2] = Orbit.candidates_image(table, phi0, th0);
            out = [];
            for ii = 1:length(k_vec)
                if (table.is_in_arc(k_vec(ii), phi1(ii)))
                    image = mod(phi1(ii), 2*pi);
                    if (abs(image - phi0) > eps)
                        out = image;
                        return;
                    end
                end
                if (table.is_in_arc(k_vec(ii), phi2(ii)))
                    image = mod(phi2(ii), 2*pi);
                    if (abs(image - phi0) > eps)
                        out = image;
                        return;
                    end
                end
            end
            phi0
            k_vec
            phi1
            phi2
            error("No potential image is the real image");
        end

        function a = arccos_argument(table, j, k, phi, th)
            rj = table.r(j);
            Oj_x = table.O(j,1);
            Oj_y = table.O(j,2);

            rk = table.r(k);
            Ok_x = table.O(k,1);
            Ok_y = table.O(k,2);
            
            aux1 = rj/rk*cos(th);
            aux2 = ((Ok_x-Oj_x)*cos(th+phi) + (Ok_y-Oj_y)*sin(th+phi))/rk;
            a =  aux1 - aux2;
        end

        function [k, phi1, phi2] = candidates_image(table, phi, th)
            % IN: table, phi, theta
            % OUT: array: (index of circle, phi1, phi2) of ray crossing
            i = table.determine_arc(phi);

            % Case we remain in the arc:
            eps = 1e-12;
            if (phi >= table.a(i) &&  phi + 2*th < table.b(i))
                k = i;
                arg = Orbit.arccos_argument(table, i, i, phi, th);
                phi1 = phi + th + acos(arg);
                phi2 = phi + th - acos(arg);
                return
            end
            if (phi >= table.a(i) &&  (phi + 2*th-table.b(i)) <= eps)
                k = mod(i, table.k) + 1;
                arg = Orbit.arccos_argument(table, i, i, phi, th);
                phi1 = phi + th + acos(arg);
                phi2 = phi + th - acos(arg);
                return
            end


            k = [];  phi1 = [];  phi2 = [];
            for jj = 1:table.k
                arg = Orbit.arccos_argument(table, i, jj, phi, th);
                if (arg > -1 && arg < 1)
                    k = [k; jj];
                    phi1 = [phi1; phi + th + acos(arg)];
                    phi2 = [phi2; phi + th - acos(arg)];
                end
            end
            if (isempty(k))
                error("No intersection of ray with circles");
            end
        end

        function obj = find_generic_sliding(table, N, phi0, th_min, th_max, max_tries)
            % finds a th0 associated to this phi0 such that the
            % corresponding orbit of N iterates is generic and sliding
            tries = 0;
            while tries < max_tries
                th0 = th_min + (th_max - th_min) * rand();
                obj = Orbit(table, phi0, th0, N);
                if ( obj.is_sliding(table) == 1 && obj.is_generic(table) == 1)
                    return;
                end
                tries = tries + 1;
            end

            error("Could not find generic sliding orbit in max_tries");
        end

        function obj = random_fundamental_generic_sliding(table, N, th_min, max_tries)
            tries = 0;
            D1 = FundamentalDomain(table, 1);
            th_max = D1.vertices(2,2);
            while tries < max_tries
                th0 =  th_min + (th_max - th_min) * rand();
                phi0 = 2 * th0 * rand();
                obj = Orbit(table, phi0, th0, N);
                if ( obj.is_sliding(table) == 1 && obj.is_generic(table) == 1)
                    return;
                end
                tries = tries + 1;
            end

            error("Could not find generic sliding orbit in max_tries");

        end

    end
end

