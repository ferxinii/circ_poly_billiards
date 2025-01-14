classdef Sequence
    %SEQUENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N  % Number of symbols
        k
        q = []  % each row represents a symbol (which has k entries)
    end

    methods 
        function plot_sequence(obj, table)
            if (obj.N ~= 4)
                error("Can only plot sequences with length 4!");
            end
            figure; 
            
            kk = 1;
            for ii = 1:4
                for jj = 1:3
                    subplot(4,4,kk);
                    title("(" + string(jj) + "," + string(obj.q(ii,jj)) + ") to (" + string(jj+1) + "," + string(obj.q(ii,jj+1)) + ")")
                    Q = FundamentalQuadrilateral(table, jj, obj.q(ii,jj));
                    Q2 = FundamentalQuadrilateral(table, jj+1, obj.q(ii, jj+1));
                    Q2.add_to_plot_lower;
                    Q2.add_to_plot_upper;
                    Q.plot_image_upper_lower(1000, table, obj.q(ii,jj));
                    kk = kk + 1; set(gca, "FontSize", 11);
                    grid on; xticklabels([]); yticklabels([]);
                end
                if (ii ~= 4)
                    subplot(4,4,kk);
                    title("(4," + string(obj.q(ii,4)) + ") to (1," + string(obj.q(ii+1,1)) + ")")
                    Q = FundamentalQuadrilateral(table, 4, obj.q(ii,4));
                    Q2 = FundamentalQuadrilateral(table, 1, obj.q(ii+1, 1));
                    Q2.add_to_plot_lower;
                    Q2.add_to_plot_upper;
                    Q.plot_image_upper_lower(1000, table, obj.q(ii,4));
                    kk = kk + 1; set(gca, "FontSize", 11);
                    grid on; xticklabels([]); yticklabels([]);
                end
            end
            sgtitle(["Properties of a finite admissible sequence (N = 4)", "Note how f stretches Q_{(j,n)} to Q_{(j+1,n')} along vertical paths"], "FontSize", 13, "FontWeight", "Bold");
        end
    end

    methods (Static)
        function obj = random_admissible(table, N, max_tries, nmin, nmax)
            obj = Sequence;
            obj.N = N;
            obj.k = table.k;
            obj.q = zeros(obj.N, obj.k);
            
            tries = 0;
            while tries < max_tries
                n = nmin; %randi(nmax) ;
                Q = FundamentalQuadrilateral(table, 1, n);
                obj.q(1,1) = n;
                for ii = 1:(N-1)
                    err = 0;
                    for jj = 2:table.k
                        try
                            n2 = Q.find_rand_n2_lemma18(table, nmax);
                            obj.q(ii,jj) = n2;
                            Q = FundamentalQuadrilateral(table, jj, n2);
                        catch
                            err = 1;
                            break;
                        end
                    end
                    if err == 1 
                        break;
                    end
                    try
                        n2 = Q.find_rand_n2_lemma18(table, nmax);
                        obj.q(ii+1,1) = n2;
                        Q = FundamentalQuadrilateral(table, 1, n2);
                    catch
                        err = 1;
                        break;
                    end
                end
                if err == 0
                    for jj = 2:table.k
                        try
                            n2 = Q.find_rand_n2_lemma18(table, nmax);
                            obj.q(end,jj) = n2;
                            Q = FundamentalQuadrilateral(table, jj, n2);
                        catch
                            err = 1;
                            break;
                        end
                    end
                    if err == 0
                        return
                    end
                end
                tries = tries + 1;
            end
            error("Could not find random admisible sequence in max_tries");
        end

        function bool = is_admissible(table, q)
            bool = 1;
            for ii = 1:(q.N-1)
                Q = FundamentalQuadrilateral(table, 1, q.q(ii,1));
                for jj = 2:q.k
                    if (~Q.pair_n_n2_satisfies_lemma18(table, q.q(ii,jj)))
                        bool = 0;
                        return;
                    end
                    Q = FundamentalQuadrilateral(table, jj, q.q(ii,jj));
                end
                if (~Q.pair_n_n2_satisfies_lemma18(table, q.q(ii+1,1)))
                    bool = 0;
                    return;
                end
            end
            
            % last row
            Q = FundamentalQuadrilateral(table, 1, q.q(end,1));
            for jj = 2:q.k
                if (~Q.pair_n_n2_satisfies_lemma18(table, q.q(end,jj)))
                    bool = 0;
                    return;
                end
                Q = FundamentalQuadrilateral(table, jj, q.q(end,jj));
            end
        end
    end
end

