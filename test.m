clear all
set(0, 'defaultFigureRenderer', 'painters')

%{  
TODO
- Proposition 29
- 
%}

%% Plot different tables
table = BilliardTable;
table.set_moss_egg(1);
table.new_plot_detailed;
title("Moss's egg", "FontSize", 13); 
subtitle("r=1", "FontSize", 13);
set(gca, "FontSize", 12);
legend();
%saveas(gcf, "../project_latex/figures/table_moss_egg.eps", "epsc");

table1 = BilliardTable;
table1.set_pseudo_ellipse(1, 2, pi/3);
table1.new_plot_detailed();
title("Pseudo-ellipse", "FontSize", 13); 
subtitle("r=1, R=2, a=\pi/3", "FontSize", 13);
set(gca, "FontSize", 12);
legend();
%saveas(gcf, "../project_latex/figures/table_pseudoellipse_1.eps", "epsc");

table2 = BilliardTable;
table2.set_pseudo_ellipse(1, 3, pi*3/4);
table2.new_plot_detailed();
title("Pseudo-ellipse", "FontSize", 13); 
subtitle("r=1, R=3, a=3\pi/4", "FontSize", 13);
set(gca, "FontSize", 12);
legend();
%saveas(gcf, "../project_latex/figures/table_pseudoellipse_2.eps", "epsc");

table3 = BilliardTable;
table3.set_table_6([1,-1], [-1,-1], [0, 1], 1);
table3.new_plot_detailed();
title("circular 6-gon", "FontSize", 13); 
subtitle("A=(1,-1), B=(-1,-1), C=(0,1), r=1", "FontSize", 13);
set(gca, "FontSize", 12);
legend();
%saveas(gcf, "../project_latex/figures/table_6gon.eps", "epsc");


%% Singularity set

L = ExtendedSingularitySet(table);

D1 = FundamentalDomain(table, 1);
D1.add_to_plot();

D1 = FundamentalDomain(table, 2);
D1.add_to_plot();

D1 = FundamentalDomain(table, 3);
D1.add_to_plot();

D1 = FundamentalDomain(table, 4);
D1.add_to_plot();

L.new_plot();

grid on;
title("Extended singularity set and fundamental domains in Moss's egg.", "FontSize", 13);
legend("D_1", "D_2", "D_3", "D_4");
set(gca, "FontSize", 12);
%saveas(gcf, "../project_latex/figures/fundamental_domains.eps", "epsc");

%%
figure;
L2 = ExtendedSingularitySet(table3);

D1 = FundamentalDomain(table3, 1);
D1.add_to_plot();

D1 = FundamentalDomain(table3, 2);
D1.add_to_plot();

D1 = FundamentalDomain(table3, 3);
D1.add_to_plot();

D1 = FundamentalDomain(table3, 4);
D1.add_to_plot();

D1 = FundamentalDomain(table3, 5);
D1.add_to_plot();

D1 = FundamentalDomain(table3, 6);
D1.add_to_plot();

L2.new_plot();

grid on;
title("Extended singularity set and fundamental domains in a 6-gon.", "FontSize", 13);
legend("D_1", "D_2", "D_3", "D_4", "D_5", "D_6");
set(gca, "FontSize", 12);
%saveas(gcf, "../project_latex/figures/fundamental_domains_6gon.eps", "epsc");


%{
%% Images of FundamentalDomain, to see Lemma 11


D1 = FundamentalDomain(table, 2);
D1.plot_sweep(table, 50, -1);
%%
D1.plot_sweep(table, 50, 1);


%%

Li = SingularitySegment;
Li = Li.j_s_neg(table, 2, 1/2);
Li.add_to_plot; hold on;
Li.plot_image(table, 1, 1000);

Li = Li.j_s_pos(table, 2, 1/2);
Li.add_to_plot;

%}

%% Example computation of the billiard map

figure;
table.add_to_plot_extended;
table.add_to_plot;
xlim("padded"); ylim("padded"); axis equal;
set(gca, "FontSize", 12);

phi = 1.5;
th = 2;

[k_vec, phi1, phi2] = Orbit.candidates_image(table, phi, th);
for ii = 1:size(k_vec,1)
    k = k_vec(ii,1);
    x = table.O(k,1) + table.r(k)*cos(phi1(ii));
    y = table.O(k,2) + table.r(k)*sin(phi1(ii));
    plot(x, y, 'o', 'LineWidth', 2, 'Color', 'blue');
    x = table.O(k,1) + table.r(k)*cos(phi2(ii));
    y = table.O(k,2) + table.r(k)*sin(phi2(ii));
    plot(x, y, 'o', 'LineWidth', 2, 'Color', 'blue');
end
k = table.determine_arc(phi);
x = table.O(k,1) + table.r(k)*cos(phi);
y = table.O(k,2) + table.r(k)*sin(phi);

plot(x, y, 'o', 'LineWidth', 2, 'Color', 'black');
xlabel("x");
ylabel("y");
legend("", "", "", "", "", "", "", "", "", "", "intersections", "initial point");

Orbit.map_phi(table, phi, th)
grid on; set(gca, "FontSize", 12);
title("Example of billiard map on Moss's egg. Initial \phi and possible images.", "FontSize", 13);
%saveas(gcf, "../project_latex/figures/map_phi.eps", "epsc");

%% Example of billiard map

phi0 = 1.2;
th0 = 0.4;

N = 100;
figure;

orbit = Orbit(table3, phi0, th0, N);
table3.add_to_plot; xlim("Padded"); ylim("Padded"); grid on;
xlabel("x"); ylabel("y");
orbit.plot_cartesian(table3); set(gca, "FontSize", 12);
axis equal; title("100 iterates of the billiard map in a 6-gon, IC = (1.2, 0.4)", "FontSize",13);
%saveas(gcf, "../project_latex/figures/map_example_1.eps", "epsc");

figure; title("100 iterates of the billiard map in a 6-gon, IC = (1.2, 0.4)", "FontSize",13);
L2.new_plot(); grid on; set(gca, "FontSize", 12);
orbit.plot_phasespace; 
%saveas(gcf, "../project_latex/figures/map_example_2.eps", "epsc");


%% Generic sliding trajectories

N = 1000;
phi0 = 0.2;
th_min = 0.2;
th_max = 1;
max_tries = 10000;

gs = Orbit.find_generic_sliding(table, N, phi0, th_min, th_max, max_tries);

sgtitle("1000 iterates of a random point of S_0^{1000} in Moss's egg", "FontSize", 13, "FontWeight", "Bold"); 

subplot(1,5,1:3)
title("Phase space");
D1 = FundamentalDomain(table, 1);
D1.add_to_plot();
D1 = FundamentalDomain(table, 2);
D1.add_to_plot();
D1 = FundamentalDomain(table, 3);
D1.add_to_plot();
D1 = FundamentalDomain(table, 4);
D1.add_to_plot();
L.add_to_plot();
gs.plot_phasespace("#0072BD");

set(gca, "FontSize", 12);
xlim([0, 2*pi]); ylim([0, 0.8])
xlabel("$\varphi$", "Interpreter", "latex");
ylabel("$\theta$", "Interpreter", "latex");
grid on;

subplot(1,5,4:5)
table.add_to_plot();
gs.plot_cartesian(table);
title("Cartesian space");
xlim("Padded"); ylim("Padded");
xlabel("x"); ylabel("y");
set(gca,'YAxisLocation','right')
axis equal; grid on;
set(gca, "FontSize", 12);

%saveas(gcf, "../project_latex/figures/generic_orbit.eps", "epsc");

%% Fundamental quadrilateral 
n = 2;
j = 1;
Q = FundamentalQuadrilateral(table, j, n);
Q.add_to_plot;

xl = xlim; yl = ylim;

Li = SingularitySegment;
Li = Li.j(table, j);
Li.add_to_plot;
Li = Li.j_s_pos(table, j, 1);
Li.add_to_plot;
Li = Li.j_s_neg(table, j+1, n);
Li.add_to_plot;
Li = Li.j_s_neg(table, j+1, n-1);
Li.add_to_plot;

legend("Q_{1,2}", "L_1", "L_1^1", "L_2^{-2}", "L_2^{-1}")

xlim([0, xl(2)]);
ylim(yl);

grid on; set(gca, "FontSize", 12);
title("(1,2)-fundamental quadrilateral", "FontSize", 13);
xlabel("$\varphi$", "Interpreter", "latex");
ylabel("$\theta$", "Interpreter", "latex");
%saveas(gcf, "../project_latex/figures/fundamental_quad_1_2.eps", "epsc");

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
grid on; set(gca, "FontSize", 12);
title("(1,2)-fundamental quadrilateral", "FontSize", 13);
xlabel("$\varphi$", "Interpreter", "latex");
ylabel("$\theta$", "Interpreter", "latex");

Li = SingularitySegment;
Li = Li.j(table, j);
Li.add_to_plot;
Li = Li.j_s_pos(table, j, 1);
Li.add_to_plot;
Li = Li.j_s_neg(table, j+1, n);
Li.add_to_plot;
Li = Li.j_s_neg(table, j+1, n-1/2);
Li.add_to_plot;
Li = Li.j_s_neg(table, j+1, n-1);
Li.add_to_plot;

legend("Q_{1,2}^+", "Q_{1,2}^-", "L_1", "L_1^1", "L_2^{-2}",  "L_2^{-1.5}", "L_2^{-1}")

xlim([0, xl(2)]);
ylim(yl);

%saveas(gcf, "../project_latex/figures/fundamental_quad_1_3.eps", "epsc");


%% Images of the fundamental quadrilateral

res = 10000;
n = 2;
j = 1;
Q = FundamentalQuadrilateral(table, j, n);

figure; 
title("f( Q_{1,2} )"); grid on; set(gca, "FontSize", 12);

Q.add_to_plot;
L.new_plot();
Q.plot_image(res, table, 1);
legend("Q_{1,2}"); set(gca, "FontSize", 12);
%saveas(gcf, "../project_latex/figures/fq_image_1.eps", "epsc");

figure;
title("f^2( Q_{1,2} )"); grid on; set(gca, "FontSize", 12);
Q.add_to_plot;
L.new_plot();
Q.plot_image(res, table, 2);
legend("Q_{1,2}"); set(gca, "FontSize", 12);
%saveas(gcf, "../project_latex/figures/fq_image_2.eps", "epsc");

figure;
title("f^3( Q_{1,2} )"); grid on; set(gca, "FontSize", 12);
Q.add_to_plot;
L.new_plot();
Q.plot_image(res, table, 3);
legend("Q_{1,2}"); set(gca, "FontSize", 12);
%saveas(gcf, "../project_latex/figures/fq_image_3.eps", "epsc");

figure;
title("f^4( Q_{1,2} )"); grid on; set(gca, "FontSize", 12);
Q.add_to_plot;
L.new_plot();
Q.plot_image(5*res, table, 4);
legend("Q_{1,2}"); set(gca, "FontSize", 12);
%saveas(gcf, "../project_latex/figures/fq_image_4.eps", "epsc");

%% Even more images of the fundamental quadrilateral

res = 10000;
n = 2;
j = 2;
Q = FundamentalQuadrilateral(table3, j, n);

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L2.new_plot; 
Q.plot_image_upper_lower(res, table3, 5);
set(gca, "FontSize", 12);
title("5 iterates of Q_{2,2} in the 6-gon", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq2_1.eps", "epsc");

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L2.new_plot; 
Q.plot_image_upper_lower(res, table3, 10);
set(gca, "FontSize", 12);
title("10 iterates of Q_{2,2} in the 6-gon", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq2_2.eps", "epsc");

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L2.new_plot; 
Q.plot_image_upper_lower(res, table3, 20);
set(gca, "FontSize", 12);
title("20 iterates of Q_{2,2} in the 6-gon", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq2_3.eps", "epsc");


figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L2.new_plot; 
Q.plot_image_upper_lower(res, table3, 30);
set(gca, "FontSize", 12);
title("30 iterates of Q_{2,2} in the 6-gon", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq2_4.eps", "epsc");


%%
res = 10000;
n = 2;
j = 1;
Q = FundamentalQuadrilateral(table, j, n);

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L.new_plot; 
Q.plot_image_upper_lower(res, table, 5);
set(gca, "FontSize", 12);
title("5 iterates of Q_{1,2} in Moss's egg", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq3_1.eps", "epsc");

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L.new_plot; 
Q.plot_image_upper_lower(res, table, 10);
set(gca, "FontSize", 12);
title("10 iterates of Q_{1,2} in Moss's egg", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq3_2.eps", "epsc");

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L.new_plot; 
Q.plot_image_upper_lower(res, table, 20);
set(gca, "FontSize", 12);
title("20 iterates of Q_{1,2} in Moss's egg", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq3_3.eps", "epsc");

figure;
Q.add_to_plot_upper;
Q.add_to_plot_lower;
L.new_plot; 
Q.plot_image_upper_lower(res, table, 30);
set(gca, "FontSize", 12);
title("30 iterates of Q_{1,2} in Moss's egg", "FontSize", 13);
saveas(gcf, "../project_latex/figures/fq3_4.eps", "epsc");

%% Exploring Lemma 17

L.new_plot();
n = 4;
j = 1;
Q = FundamentalQuadrilateral(table, j, n);
Q.add_to_plot;
Q.plot_domain_upper_lower(1000);
Q.plot_image_upper_lower(10000, table, n);

figure;
L.new_plot();
n = 100;
j = 1;
Q = FundamentalQuadrilateral(table, j, n);
Q.add_to_plot;
Q.plot_domain_upper_lower(1000);
Q.plot_image_upper_lower(10000, table, n);

% However, now muj > 1 !!
figure;
L.new_plot();
n = 4;
j = 2;
Q = FundamentalQuadrilateral(table, j, n);
Q.add_to_plot;
Q.plot_domain_upper_lower(1000);
Q.plot_image_upper_lower(10000, table, n);

%% Lemma 17

j = 3;
n = 5;
s = 0;
Q = FundamentalQuadrilateral(table, j, n);
Q.add_to_plot;
L.add_to_plot;
L1 = SingularitySegment;
L1 = L1.j_s_neg(table, j+1, n-1/2);
L1.add_to_plot;
L.add_to_plot;
Lemma17(table, j, n)

figure;
j = 2;
n = 5;
s = 0;
Q = FundamentalQuadrilateral(table, j, n);
Q.add_to_plot;
L.add_to_plot;
L1 = SingularitySegment;
L1 = L1.j_s_neg(table, j+1, n-1/2);
L1.add_to_plot;
L.add_to_plot;
Lemma17(table, j, n)

%% Lemma 18
% OBS: Lemma talks about stretching vertical paths into vertical paths,
% but not that the orientations bust be the same!
nmax = 100;

j = 1;
n = 30;
Q = FundamentalQuadrilateral(table, j, n);
n2 = Q.find_min_n2_lemma18(table, nmax);
Q2 = FundamentalQuadrilateral(table, mod(j, table.k) + 1, n2);
Q2.add_to_plot_lower;
Q2.add_to_plot_upper;
Q.plot_image_upper_lower(5000, table, n);
L.add_to_plot;
grid on;
xlabel("$\varphi$", "Interpreter", "latex");
ylabel("$\theta$", "Interpreter", "latex");
set(gca, "FontSize", 12);
title("f^{30}(Q_{1,30}) on top of Q_{2," + string(n2) + "} in Moss's egg.", "FontSize", 13);
saveas(gcf, "../project_latex/figures/stretching_along_paths_1.eps", "epsc");
mu1 = sqrt(table.r(1) / table.r(2))


figure;
j = 2;
n = 7;
Q = FundamentalQuadrilateral(table, j, n);
n2 = Q.find_min_n2_lemma18(table, nmax);
Q2 = FundamentalQuadrilateral(table, mod(j, table.k) + 1, n2);
Q2.add_to_plot_lower;
Q2.add_to_plot_upper;
Q.plot_image_upper_lower(5000, table, n);
L.add_to_plot;
grid on; 
xlabel("$\varphi$", "Interpreter", "latex");
ylabel("$\theta$", "Interpreter", "latex");
set(gca, "FontSize", 12);
title("f^{7}(Q_{2,7}) on top of Q_{3," + string(n2) + "} in Moss's egg.", "FontSize", 13);
saveas(gcf, "../project_latex/figures/stretching_along_paths_2.eps", "epsc");
mu2 = sqrt(table.r(2) / table.r(3))


%% Admissible sequences

nmin = 35;
nmax = 100;
max_tries = 1000;
N = 4;

%q = Sequence.random_admissible(table, N, max_tries, nmin, nmax);
Sequence.is_admissible(table, q)

q.plot_sequence(table)

%saveas(gcf, "../project_latex/figures/admissible_sequence.eps", "epsc");

%% Return map

N = 10000;
max_tries = 1000;
th_min = 0.0001;
orbit = Orbit.random_fundamental_generic_sliding(table, N, th_min, max_tries);
F = orbit.extract_return_map(table)

F.plot_phasespace()


