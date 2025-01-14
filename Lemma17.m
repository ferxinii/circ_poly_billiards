classdef Lemma17
    %LEMMA17 Summary of this class goes here
    
    methods 
        function obj = Lemma17(table, j, n)
            if (n < table.chi_min(j))
                table.chi_min(j)
                error("n must be >= than chi_min(j) (to compare with lemma 17)!");
            end
            disp("xhi(j) = " + num2str(table.chi_min(j), "%d"));

            res = 1000;
            [nu_a, w_a] = Lemma17.nuw_j_n(table, j, n);
            [nu_b0, w_b0] = Lemma17.nuw_j_n_s(res, table, j, n, 0);
            [nu_b1, w_b1] = Lemma17.nuw_j_n_s(res, table, j, n, 1);
            [nu_b12, w_b12] = Lemma17.nuw_j_n_s(res, table, j, n, 1/2);
            
            dj = table.b(j) - table.a(j);
            rj = table.r(j);
            rnext = table.r(mod(j, table.k) + 1);
            muj = sqrt(rj / rnext);
            
            format long
            disp("a)");
            disp("nu(j,n) = " + num2str(nu_a, "%.16f") + "  VS  " + num2str(dj/(2*n+2), "%.16f"));
            disp("w(j,n) = " + num2str(w_a, "%.16f") + "  VS  " + num2str(dj/(2*n-2), "%.16f"));
            disp("b)");
            disp("nu(j,n,0) = " + num2str(nu_b0, "%.16f") + "  VS  " + num2str(dj/(2*n+2), "%.16f"));
            disp("w(j,n,0) = " + num2str(w_b0, "%.16f") + "  VS  " + num2str(dj/(2*n), "%.16f"));
            disp("nu(j,n,1) = " + num2str(nu_b1, "%.16f") + "  VS  " + num2str(dj/(2*n), "%.16f"));
            disp("w(j,n,1) = " + num2str(w_b1, "%.16f") + "  VS  " + num2str(dj/(2*n-2), "%.16f"));
            if (muj < 1)
                disp("muj = " + num2str(muj) + " < 1:  w(j,n,1/2) = " + num2str(w_b12, "%.16f") + "  VS  " + num2str(muj*dj/(2*n-1), "%.16f"));
            elseif (muj > 1)
                disp("muj = " + num2str(muj) + " > 1:  nu(j,n,1/2) = " + num2str(nu_b12, "%.16f") + "  VS  " + num2str(muj*dj/(2*n+1), "%.16f"));
            else 
                disp("muj = 1...");
            end
        end
    end

    methods (Static)
        function [nu, w] = nuw_j_n(table, j, n)
            Q = FundamentalQuadrilateral(table, j, n);
            nu = min(Q.vertices(:,2));
            w = max(Q.vertices(:,2));
        end

        function [nu, w] = nuw_j_n_s(resolution, table, j, n, s)
            % Lemma 17
            if (n < table.chi_min(j))
                table.chi_min(j)
                error("n must be >= than chi_min(j) (to compare with lemma 17)!");
            end
            if (s < 0 || s > 1)
                error("domain is empty!")
            end
            [xi, yi] = SingularitySegment.intersection_segments(table, j, 0, n-s);
            [xf, yf] = SingularitySegment.intersection_segments(table, j, 1, n-s);
            
            domain = linspace_vec([xi, yi], [xf, yf], resolution)';
            image = zeros(size(domain));
            for ii = 1:length(domain)
                orbit = Orbit(table, domain(ii,1), domain(ii,2), n);
                image(ii,:) = orbit.iter(end,:);
            end
            nu = min(image(:, 2));
            w = max(image(:, 2));

            plot([domain(:,1); image(:,1)], [domain(:,2); image(:,2)], "x"); 
        end
    end
end

function arr = linspace_vec(vi, vf, N)
    arr = bsxfun(@plus,((vf(:)-vi(:))./(N-1))*[0:N-1],vi(:));
end