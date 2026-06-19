%% spiral_quadratic_trapping.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: Quadratic trapping regions for a spiral attractor
%
% This script computes quadratic trapping regions for the planar system
%
%     x' = x(1-r^2) - omega*y,
%     y' = y(1-r^2) + omega*x,
%
% where r^2 = x^2 + y^2. The origin is unstable and the unit circle is an
% attracting periodic orbit.
%
% For several values of lambda, the script searches for a positive definite
% matrix P such that the set
%
%     { z : z' P z <= 1 }
%
% is trapping. The SOS condition is
%
%     lambda * grad(V) * f <= -V,
%
% where V(z) = z' P z - 1.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     spiral_sos_quadratic_trapping_boundaries.pdf
%
% in the current directory.
%
% Requirements:
%   - MATLAB
%   - YALMIP
%   - MOSEK
%
% -------------------------------------------------------------------------

clear; clc; close all;
yalmip('clear');

%% User options

omega = 2.0;

lambda_vals = [0.10 0.25 0.50 1.00 2.00];

eps_pd = 1e-8;

solver_name = 'mosek';
verbose_solver = 0;

saveFigure = true;
figureFile = 'spiral_sos_quadratic_trapping_boundaries.pdf';

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',1);

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.orange = [230, 120, 20]/255;
S.purple = [170, 90, 160]/255;
S.cyan   = [0, 170, 200]/255;
S.black  = [0, 0, 0];

region_colors = {S.blue, S.red, S.orange, S.purple, S.cyan};

%% Symbolic variables and vector field

sdpvar x y
z = [x; y];

r2 = x^2 + y^2;

f = [
    x*(1-r2) - omega*y;
    y*(1-r2) + omega*x
];

%% SOS computation

fprintf('\n============================================================\n');
fprintf('Quadratic trapping regions for a spiral attractor\n');
fprintf('============================================================\n\n');

fprintf('omega = %.6g\n\n', omega);

results = struct();

for j = 1:length(lambda_vals)

    lambda = lambda_vals(j);

    sdpvar p11 p12 p22
    P = [p11 p12; p12 p22];

    V = z.'*P*z - 1;
    dVf = jacobian(V,z)*f;

    G = -V - lambda*dVf;

    cons = [
        P >= eps_pd*eye(2), ...
        sos(G)
        ];

    % Feasibility is enough for the example. The normalization is implicit
    % in the level set z'Pz <= 1.
    sol = solvesos(cons, [], opts, [p11; p12; p22]);

    results(j).lambda = lambda;
    results(j).feasible = (sol.problem == 0);
    results(j).info = sol.info;

    if results(j).feasible

        results(j).P = value(P);
        results(j).area = pi/sqrt(det(results(j).P));

        fprintf('lambda = %.3f feasible\n', lambda);
        fprintf('P = \n');
        disp(results(j).P);
        fprintf('area = %.8f\n\n', results(j).area);

    else

        results(j).P = nan(2);
        results(j).area = nan;

        fprintf('lambda = %.3f infeasible: %s\n\n', lambda, sol.info);
    end
end

%% Plotting window and vector field

xmin = -2.1; xmax = 2.1;
ymin = -2.1; ymax = 2.1;

nx = 17;
ny = 17;

[Xg,Yg] = meshgrid(linspace(xmin,xmax,nx), linspace(ymin,ymax,ny));

Rg2 = Xg.^2 + Yg.^2;

Ug = Xg.*(1 - Rg2) - omega*Yg;
Vg = Yg.*(1 - Rg2) + omega*Xg;

Spd = sqrt(Ug.^2 + Vg.^2);
Ug = Ug ./ max(Spd,1e-12);
Vg = Vg ./ max(Spd,1e-12);

%% Sample trajectories

ICs = [
     0.20   0.00
     0.55   0.25
    -0.45   0.55
     1.65   0.10
    -1.55  -0.55
     0.15  -1.70
     1.25  -1.35
];

tspan = [0 18];
odeopts = odeset('RelTol',1e-10,'AbsTol',1e-12);

theta = linspace(0,2*pi,1000);

%% Plot

fig = figure;
set(fig,'Color','w','Units','centimeters','Position',[2 2 11 11]);

hold on
box on

quiver(Xg,Yg,Ug,Vg,0.30, ...
    'Color', 0.86*[1 1 1], ...
    'LineWidth', 0.6, ...
    'MaxHeadSize', 0.55);

for k = 1:size(ICs,1)

    z0 = ICs(k,:)';

    [~,Ztraj] = ode45(@(t,z) spiral_rhs(t,z,omega), ...
        tspan, z0, odeopts);

    plot(Ztraj(:,1), Ztraj(:,2), '-', ...
        'Color', 0.65*S.blue + 0.35*S.black, ...
        'LineWidth', 3);

    plot(Ztraj(1,1), Ztraj(1,2), 'o', ...
        'MarkerSize', 6.5, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', S.black, ...
        'LineWidth', 1);
end

h_unit = plot(cos(theta), sin(theta), '-', ...
    'Color', S.black, ...
    'LineWidth', 4);

h_trap = gobjects(length(lambda_vals),1);
leg_entries = cell(length(lambda_vals),1);

for j = 1:length(lambda_vals)

    if ~results(j).feasible
        continue;
    end

    P = results(j).P;
    [xe, ye] = ellipse_from_P(P, theta);

    plot(xe, ye, '--', ...
        'Color', region_colors{j}, ...
        'LineWidth', 4);

    h_trap(j) = plot(nan,nan,'--', ...
        'Color', region_colors{j}, ...
        'LineWidth', 4);

    leg_entries{j} = sprintf('$\\lambda=%.2f$', results(j).lambda);
end

h_eq = plot(0,0,'o', ...
    'MarkerSize', 20, ...
    'MarkerFaceColor', S.green, ...
    'MarkerEdgeColor', S.black, ...
    'LineWidth', 1);

axis equal
axis([xmin xmax ymin ymax]);

set(gca,'FontSize',24)
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold');

valid = isgraphics(h_trap);
h_leg = [h_trap(valid); h_unit; h_eq];
leg_labels = [leg_entries(valid); {'unit circle'; 'equilibrium'}];

legend(h_leg, leg_labels, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'FontSize',18, ...
    'Box','on');

if saveFigure
    save_pdf_figure(fig,figureFile);
end

fprintf('Finished spiral trapping-region example.\n');

%% Local functions

function dzdt = spiral_rhs(~,z,omega)

    x = z(1);
    y = z(2);

    r2 = x^2 + y^2;

    dzdt = [
        x*(1-r2) - omega*y;
        y*(1-r2) + omega*x
    ];
end

function [xe, ye] = ellipse_from_P(P, theta)
% Boundary of {z'Pz <= 1}.

    [Q,D] = eig(P);

    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));

    E = Q*[a*cos(theta); b*sin(theta)];

    xe = E(1,:);
    ye = E(2,:);
end

function save_pdf_figure(fig,fileName)

    set(fig, 'PaperUnits','centimeters');
    set(fig, 'Units','centimeters');

    pos = get(fig,'Position');

    set(fig, 'PaperSize', [pos(3) pos(4)]);
    set(fig, 'PaperPositionMode', 'manual');
    set(fig, 'PaperPosition', [0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n', fileName);
end