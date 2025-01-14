classdef ExtendedSingularitySet
    % The extended singularity set associated to a particular table
    
    properties
        L = [SingularitySegment]  % Array of SignularitySegment
    end
    
    methods
        function obj = ExtendedSingularitySet(table)
            %CONSTRUCTOR
            L_aux = SingularitySegment;
            obj.L(table.k * 3, 1) = L_aux;  % initialize array of objects
            for ii = 1:table.k
                obj.L((ii-1)*3 + 1) = L_aux.j(table, ii);
                obj.L((ii-1)*3 + 2) = L_aux.j_s_pos(table, ii, 1/2);
                obj.L((ii-1)*3 + 3) = L_aux.j_s_pos(table, ii, 1);
            end
        end
        
        function new_plot(obj)
             colors = {
                [0.2, 0.2, 0.2];  % Dark Gray
                [0.5, 0.5, 0.5];  % Gray
                [0.7, 0.7, 0.7];  % Light Gray
            };
            hold on;
            set(gca, "FontSize", 12);
            xlabel("$\varphi$", "Interpreter", "latex");
            ylabel("$\theta$", "Interpreter", "latex");
            xlim([0, 2*pi]); ylim([0, pi]);
            for ii = 1:length(obj.L)
                obj.L(ii).add_to_plot(colors{mod(ii-1,3) + 1});
            end
        end

        function add_to_plot(obj)
            colors = {
                [0.2, 0.2, 0.2];  % Dark Gray
                [0.5, 0.5, 0.5];  % Gray
                [0.7, 0.7, 0.7];  % Light Gray
            };
            xl = xlim(); yl = ylim();
            hold on;
            for ii = 1:length(obj.L)
                obj.L(ii).add_to_plot(colors{mod(ii-1,3) + 1});
            end
            xlim(xl); ylim(yl); 
        end
    end
end

