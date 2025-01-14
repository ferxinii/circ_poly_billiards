classdef FundamentalQuadrilateral
    %FUNDAMENTALQUADRILATERAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        j
        n
        vertices = zeros(4,2);
        upper =  zeros(4,2);
        lower =  zeros(4,2);
    end
    
    methods
        function obj = FundamentalQuadrilateral(table, j, n)
            if (n < 2 || j < 1 || j > table.k)
                error("Bad indices for fundamental quadrilateral");
            end
            obj.j = j;
            obj.n = n;

            % By Lemma 15
            aj = table.a(j);
            dj = table.b(j) - table.a(j);

            obj.vertices(1, :) = [aj, dj/(2*n-2)];
            obj.vertices(2, :) = [aj, dj/(2*n)];
            obj.vertices(3, :) = [aj+dj/(n+1), dj/(2*n + 2)];
            obj.vertices(4, :) = [aj+dj/n, dj/(2*n)];

            obj.lower(1, :) = [aj, dj/(2*n-1)];
            obj.lower(2, :) = [aj, dj/(2*n)];
            obj.lower(3, :) = [aj+dj/(n+1), dj/(2*n + 2)];
            obj.lower(4, :) = [aj+dj/(n+1/2), dj/(2*n+1)];

            obj.upper(1, :) = [aj, dj/(2*n-2)];
            obj.upper(2, :) = [aj, dj/(2*n-1)];
            obj.upper(3, :) = [aj+dj/(n+1/2), dj/(2*n + 1)];
            obj.upper(4, :) = [aj+dj/n, dj/(2*n)];
        end
        
        function add_to_plot(obj)
            hold on;
            plot(polyshape(obj.vertices));
        end

        function add_to_plot_lower(obj)
            hold on;
            plot(polyshape(obj.lower), "FaceColor", "blue");
        end

        function add_to_plot_upper(obj)
            hold on;
            plot(polyshape(obj.upper), "FaceColor", "red");
        end

        function image = image_quadr(obj, vertices, resolution, eps, table, n_it)
            % eps adds some margin towards the interior.
            N = resolution / 4;
            
            xi = vertices(1,1) + eps;
            yi = vertices(1,2) - eps;
            xf = vertices(2,1) + eps;
            yf = vertices(2,2) + eps;
            left = linspace_vec([xi, yi], [xf, yf], N)';

            xi = vertices(2,1) + eps;
            yi = vertices(2,2) + eps;
            xf = vertices(3,1) - eps;
            yf = vertices(3,2) + eps;
            bot = linspace_vec([xi, yi], [xf, yf], N)';

            xi = vertices(3,1) - eps;
            yi = vertices(3,2) + eps;
            xf = vertices(4,1) - eps;
            yf = vertices(4,2) - eps;
            right = linspace_vec([xi, yi], [xf, yf], N)';

            xi = vertices(4,1) - eps;
            yi = vertices(4,2) - eps;
            xf = vertices(1,1) + eps;
            yf = vertices(1,2) - eps;
            top = linspace_vec([xi, yi], [xf, yf], N)';

            sides = [left; bot; right; top];
            %plot(sides(:,1), sides(:,2), "o");  % DEBUG
            image = zeros(size(sides));
            for ii = 1:length(sides)
                orbit = Orbit(table, sides(ii,1), sides(ii,2), n_it);
                image(ii,:) = orbit.iter(end,:);
            end
        end

        function plot_image(obj, resolution, table, n_it)
            eps = 1e-14;
            
            hold on;
            im = obj.image_quadr(obj.vertices, resolution, eps, table, n_it);
            plot(im(:,1), im(:,2), ".");
        end

        function plot_image_upper_lower(obj, resolution, table, n_it)
            eps = 1e-10;
            
            hold on;
            im = obj.image_quadr(obj.upper, resolution, eps, table, n_it);
            plot(im(:,1), im(:,2), ".", "Color", "red");

            im = obj.image_quadr(obj.lower, resolution, eps, table, n_it);
            plot(im(:,1), im(:,2), ".", "Color", "blue");
        end

        function plot_domain(obj, resolution)
            hold on;
            N = resolution/4;
            left = linspace_vec(obj.vertices(1,:), obj.vertices(2,:), N)';
            bot = linspace_vec(obj.vertices(2,:), obj.vertices(3,:), N)';
            right = linspace_vec(obj.vertices(3,:), obj.vertices(4,:), N)';
            top = linspace_vec(obj.vertices(4,:), obj.vertices(1,:), N)';
            sides = [left; bot; right; top];
            plot(sides(:,1), sides(:,2), ".");

        end

        function plot_domain_upper_lower(obj, resolution)
            N = resolution/4;
            left = linspace_vec(obj.upper(1,:), obj.upper(2,:), N)';
            bot = linspace_vec(obj.upper(2,:), obj.upper(3,:), N)';
            right = linspace_vec(obj.upper(3,:), obj.upper(4,:), N)';
            top = linspace_vec(obj.upper(4,:), obj.upper(1,:), N)';
            sides = [left; bot; right; top];
            plot(sides(:,1), sides(:,2), ".", "Color", "red");


            left = linspace_vec(obj.lower(1,:), obj.lower(2,:), N)';
            bot = linspace_vec(obj.lower(2,:), obj.lower(3,:), N)';
            right = linspace_vec(obj.lower(3,:), obj.lower(4,:), N)';
            top = linspace_vec(obj.lower(4,:), obj.lower(1,:), N)';
            sides = [left; bot; right; top];
            plot(sides(:,1), sides(:,2), ".", "Color", "blue");
        end

        function bool = pair_n_n2_satisfies_lemma18(obj, table, n2)
            % This cheks the conditions of lemma 18
            if (obj.n < table.chi_min(obj.j))
                bool = 0;
                return;
            end

            if (n2 < table.chi_min(mod(obj.j, table.k) + 1))
                bool = 0;
                return;
            end
            
            dj = table.b(obj.j) - table.a(obj.j);
            dnext = table.b(mod(obj.j, table.k) + 1) - table.a(mod(obj.j, table.k) + 1);
            rj = table.r(obj.j);
            rnext = table.r(mod(obj.j, table.k) + 1);
            muj = sqrt(rj / rnext);
            
            ajm = dnext / (dj * max(1, muj));
            ajp = dnext / (dj * min(1, muj));
            bjp = ajp + 1;
            bjm = ajm + 1;

            cond1 = ajm * obj.n + bjm <= n2;
            cond2 = n2 <= ajp * obj.n - bjp;
            if (cond1 && cond2)
                bool = 1;
            else
                bool = 0;
            end
        end

        function n2 = find_min_n2_lemma18(obj, table, nmax)
            for n2 = 1:nmax
                if (obj.pair_n_n2_satisfies_lemma18(table, n2) == 1)
                    return;
                end
            end
            error("could not find n2 <= nmax!, fix this function to not give error?");
        end

        function n2 = find_rand_n2_lemma18(obj, table, nmax)
            tries = 0;
            while (tries < nmax)
                n2 = randi(nmax);
                if (obj.pair_n_n2_satisfies_lemma18(table, n2) == 1)
                    return;
                end
                tries = tries + 1;
            end
            error("could not find n2 <= nmax!, fix this function to not give error?");
        end
    end
        
end

function arr = linspace_vec(vi, vf, N)
    arr = bsxfun(@plus,((vf(:)-vi(:))./(N-1))*[0:N-1],vi(:));
end


