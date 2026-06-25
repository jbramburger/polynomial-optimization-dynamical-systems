%% henon_heiles_lyapunov_exponent.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Maximal Lyapunov exponent bounds for the Henon-Heiles system
%
% This script computes SOS upper bounds on the maximal Lyapunov exponent
% for the Henon-Heiles system and produces the associated energy-set and
% extremal-orbit figures.
%
% Pressing Run can reproduce:
%
%   henon_heiles_le_bounds.csv
%   henon_heiles_energy_projection.pdf
%   henon_heiles_extremal_orbit_3d.pdf
%   henon_heiles_extremal_orbit_3d_x4.pdf
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

runSOS = true;          % Set false to regenerate only the figures.
saveFigures = true;
saveTable = true;

degrees = [2 4 6 8 10];

solver_name = 'mosek';
verbose_solver = 1;

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',1);

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.black  = [0, 0, 0];
S.light_gray = [0.78, 0.78, 0.78];
S.fill_blue = [120, 170, 220]/255;

fprintf('\n============================================================\n');
fprintf('Henon-Heiles maximal Lyapunov exponent via SOS\n');
fprintf('============================================================\n\n');

%% Part I: SOS hierarchy

if runSOS

    Bounds = nan(numel(degrees),1);
    Status = strings(numel(degrees),1);

    for k = 1:numel(degrees)

        d = degrees(k);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('SOS Lyapunov-exponent bound, degree d = %d\n', d);
        fprintf('------------------------------------------------------------\n');

        [Bval, info] = solve_henon_heiles_le_bound(d, opts);

        Bounds(k) = Bval;
        Status(k) = string(info);

        fprintf('  Bound = %.12g   [%s]\n', Bval, info);
    end

    T = table(degrees(:), Bounds, Status, ...
        'VariableNames', {'Degree','UpperBound','Status'});

    disp(T);

    if saveTable
        writetable(T,'henon_heiles_le_bounds.csv');
        fprintf('Saved table: henon_heiles_le_bounds.csv\n');
    end
end

%% Part II: energy projection and extremal periodic orbit

plot_henon_heiles_energy_and_orbits(S, saveFigures);

fprintf('\nFinished Henon-Heiles Lyapunov-exponent example.\n');

%% Local functions

function [Bval, solinfo] = solve_henon_heiles_le_bound(d, opts)

    yalmip('clear');

    sdpvar x1 x2 x3 x4
    sdpvar z1 z2 z3 z4
    sdpvar B

    x = [x1; x2; x3; x4];
    z = [z1; z2; z3; z4];
    vars = [x; z];

    H = 0.5*(x1^2 + x2^2 + x3^2 + x4^2) ...
        + x1^2*x2 - (1/3)*x2^3;

    f = [
        x3;
        x4;
       -x1 - 2*x1*x2;
       -x2 - x1^2 + x2^2
    ];

    Df = jacobian(f,x);

    Phi = z.'*Df*z;
    ell = Df*z - Phi*z;

    [V, cV]       = polynomial(vars,d);
    [sigma1, cS1] = polynomial(vars,d);
    [sigma2, cS2] = polynomial(vars,d);
    [sigma3, cS3] = polynomial(vars,d);
    [rho0, cR0]   = polynomial(vars,d);

    g1 = 1/7 - H;
    g2 = H;
    g3 = 1 - x1^2 - x2^2;

    h0 = 1 - (z1^2 + z2^2 + z3^2 + z4^2);

    P = B - Phi - jacobian(V,x)*f - jacobian(V,z)*ell;

    residual = P - sigma1*g1 - sigma2*g2 - sigma3*g3 - rho0*h0;

    constraints = [
        sos(residual), ...
        sos(sigma1), ...
        sos(sigma2), ...
        sos(sigma3)
        ];

    Gens = henon_heiles_symmetry_generators();

    tunablePolys = {V, sigma1, sigma2, sigma3, rho0};

    for q = 1:numel(tunablePolys)
        constraints = [constraints, ...
            invariant_constraints(tunablePolys{q}, vars, Gens)]; 
    end

    decvars = [B; cV; cS1; cS2; cS3; cR0];

    sol = solvesos(constraints, B, opts, decvars);

    Bval = value(B);
    solinfo = sol.info;

    if sol.problem ~= 0
        warning('Solver status at degree %d: %s', d, sol.info);
    end
end

function Gens = henon_heiles_symmetry_generators()

    Ref2 = [-1 0; 0 1];

    theta = 2*pi/3;
    Rot2 = [cos(theta) -sin(theta);
            sin(theta)  cos(theta)];

    Ref4 = blkdiag(Ref2,Ref2);
    Rot4 = blkdiag(Rot2,Rot2);

    G_reflect = blkdiag(Ref4,Ref4);
    G_rotate  = blkdiag(Rot4,Rot4);
    G_tangent_sign = blkdiag(eye(4),-eye(4));

    Gens = {G_reflect, G_rotate, G_tangent_sign};
end

function Fsym = invariant_constraints(p, vars, Gens)

    Fsym = [];

    for j = 1:numel(Gens)

        G = Gens{j};
        pG = replace(p, vars, G*vars);
        diffPoly = clean(p - pG, 1e-12);

        coeffs = coefficients(diffPoly, vars);
        Fsym = [Fsym, coeffs == 0];
    end
end

function plot_henon_heiles_energy_and_orbits(S, saveFigures)

    Emax = 1/7;

    U = @(q1,q2) 0.5*(q1.^2 + q2.^2) ...
        + q1.^2.*q2 - (1/3)*q2.^3;

    N = 800;
    q1min = -1.15; q1max = 1.15;
    q2min = -1.15; q2max = 1.15;

    [q1,q2] = meshgrid(linspace(q1min,q1max,N), ...
                       linspace(q2min,q2max,N));

    Uval = U(q1,q2);

    allowedEnergy = Uval <= Emax;

    %% Figure 1: configuration-space projection

    fig1 = figure;
    set(fig1,'Color','w','Units','centimeters','Position',[2 2 12 11]);

    hold on
    box on

    contourf(q1, q2, double(allowedEnergy), [0.5 1.5], ...
        'LineStyle','none', ...
        'FaceAlpha',0.95);

    colormap(gca, [1 1 1; S.fill_blue]);

    contour(q1, q2, Uval, [Emax Emax], ...
        'Color', S.blue, ...
        'LineWidth', 5);

    innerLevels = Emax*[0.15 0.35 0.55 0.75];

    contour(q1, q2, Uval, innerLevels, ...
        'Color', S.light_gray, ...
        'LineWidth', 5);

    th = linspace(0,2*pi,600);

    plot(cos(th), sin(th), '--', ...
        'Color', S.red, ...
        'LineWidth', 5);

    axis equal
    axis([q1min q1max q2min q2max]);

    xlabel('$x_1$','Interpreter','latex','FontSize',42,'FontWeight','bold');
    ylabel('$x_2$','Interpreter','latex','FontSize',42,'FontWeight','bold');

    set(gca,'FontSize',42,'TickLabelInterpreter','latex','Layer','top');

    if saveFigures
        save_pdf_figure(fig1,'henon_heiles_energy_projection.pdf');
    end

    %% Periodic orbit

    x0 = [
        0.562878385826716;
       -0.053847890920149;
        0;
        0
    ];

    Tper = 6.966517640959103;

    f = @(~,x) [
        x(3);
        x(4);
       -x(1) - 2*x(1)*x(2);
       -x(2) - x(1)^2 + x(2)^2
    ];

    optsODE = odeset('RelTol',1e-11,'AbsTol',1e-13);

    [~,u] = ode113(f, linspace(0,Tper,2000), x0, optsODE);

    fprintf('Initial point error after one period: %.3e\n', norm(u(end,:).'-x0));
    fprintf('Energy at x0: %.15f\n', Hval(x0));
    fprintf('Energy variation along orbit: %.3e\n', ...
        max(abs(arrayfun(@(i) Hval(u(i,:).'), 1:size(u,1)) - Hval(x0))));

    theta = 2*pi/3;

    RotP = [cos(theta) -sin(theta);
            sin(theta)  cos(theta)];

    RotM = [cos(-theta) -sin(-theta);
            sin(-theta)  cos(-theta)];

    R4P = blkdiag(RotP,RotP);
    R4M = blkdiag(RotM,RotM);

    uP = (R4P*u.').';
    uM = (R4M*u.').';

    %% Figure 2: 3D orbit in (x1,x2,x3)

    fig2 = figure;
    set(fig2,'Color','w','Units','centimeters','Position',[2 2 12 11]);

    hold on
    box on
    grid on

    plot3(u(:,1),  u(:,2),  u(:,3),  '-', 'Color', S.green, 'LineWidth', 8);
    plot3(uP(:,1), uP(:,2), uP(:,3), '-', 'Color', S.blue,  'LineWidth', 8);
    plot3(uM(:,1), uM(:,2), uM(:,3), '-', 'Color', S.red,   'LineWidth', 8);

    xlabel('$x_1$','Interpreter','latex','FontSize',42,'FontWeight','bold');
    ylabel('$x_2$','Interpreter','latex','FontSize',42,'FontWeight','bold');
    zlabel('$x_3$','Interpreter','latex','FontSize',42,'FontWeight','bold');

    set(gca,'FontSize',42,'TickLabelInterpreter','latex');
    axis tight
    view(38,24);

    if saveFigures
        save_pdf_figure(fig2,'henon_heiles_extremal_orbit_3d.pdf');
    end

    %% Figure 3: 3D orbit in (x1,x2,x4)

    fig3 = figure;
    set(fig3,'Color','w','Units','centimeters','Position',[2 2 12 11]);

    hold on
    box on
    grid on

    plot3(u(:,1),  u(:,2),  u(:,4),  '-', 'Color', S.green, 'LineWidth', 8);
    plot3(uP(:,1), uP(:,2), uP(:,4), '-', 'Color', S.blue,  'LineWidth', 8);
    plot3(uM(:,1), uM(:,2), uM(:,4), '-', 'Color', S.red,   'LineWidth', 8);

    xlabel('$x_1$','Interpreter','latex','FontSize',42,'FontWeight','bold');
    ylabel('$x_2$','Interpreter','latex','FontSize',42,'FontWeight','bold');
    zlabel('$x_4$','Interpreter','latex','FontSize',42,'FontWeight','bold');

    set(gca,'FontSize',42,'TickLabelInterpreter','latex');
    axis tight
    view(38,24);

    if saveFigures
        save_pdf_figure(fig3,'henon_heiles_extremal_orbit_3d_x4.pdf');
    end
end

function H = Hval(x)

    H = 0.5*(x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2) ...
        + x(1)^2*x(2) - (1/3)*x(2)^3;
end

function save_pdf_figure(fig,fileName)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n', fileName);
end