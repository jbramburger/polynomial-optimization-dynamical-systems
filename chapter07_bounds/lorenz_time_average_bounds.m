%% lorenz_time_average_bounds_sos.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Time-average bounds for normalized moments of the Lorenz system
%
% This script computes SOS auxiliary-function upper bounds for normalized
% moments of the Lorenz system on a compact rescaled box.
%
% Pressing Run should reproduce the numerical output and save:
%
%     lorenz_normalized_moment_bounds.csv
%     lorenz_z4_defect_xz.pdf
%     lorenz_z4_defect_xy.pdf
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

sigma = 10;
rho   = 28;
beta  = 8/3;

degrees = [4 6 8 10 12];

solver_name = 'mosek';
verbose_solver = 1;

saveFigures = true;
saveTable = true;

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',1);

%% Book plotting style

S.blue  = [0, 92, 175]/255;
S.red   = [200, 50, 50]/255;
S.black = [0, 0, 0];

%% Rescaled box K = [-1,1]^3

Kbox.mid = [0; 0; 25];
Kbox.rad = [25; 35; 25];

%% Observables

obs = {
    'z',        [0 0 1]
    'x^2',      [2 0 0]
    'xy',       [1 1 0]
    'y^2',      [0 2 0]
    'z^2',      [0 0 2]
    'x^2 z',    [2 0 1]
    'y^2 z',    [0 2 1]
    'xyz',      [1 1 1]
    'z^3',      [0 0 3]
    'x^4',      [4 0 0]
    'x^3 y',    [3 1 0]
    'x^2 y^2',  [2 2 0]
    'x^2 z^2',  [2 0 2]
    'xy^3',     [1 3 0]
    'xyz^2',    [1 1 2]
    'y^4',      [0 4 0]
    'y^2 z^2',  [0 2 2]
    'z^4',      [0 0 4]
};

%% Run SOS bounds

Bounds = nan(size(obs,1), numel(degrees));
Status = strings(size(obs,1), numel(degrees));

fprintf('\n============================================================\n');
fprintf('Lorenz time-average bounds via auxiliary functions\n');
fprintf('============================================================\n\n');

for i = 1:size(obs,1)

    name  = obs{i,1};
    expnt = obs{i,2};
    obsDeg = sum(expnt);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Observable: %s\n', name);
    fprintf('------------------------------------------------------------\n');

    for k = 1:numel(degrees)

        d = degrees(k);

        if d < obsDeg
            fprintf('  d = %d skipped: observable degree is %d.\n', d, obsDeg);
            Status(i,k) = "skipped";
            continue
        end

        degV = d - 1;

        fprintf('  Solving d = %d, deg(V) <= %d ...\n', d, degV);

        [Cval, solinfo] = solve_lorenz_moment_bound( ...
            expnt, d, degV, sigma, rho, beta, Kbox, opts);

        Bounds(i,k) = Cval;
        Status(i,k) = string(solinfo);

        fprintf('    C_%d = %.12g   [%s]\n', d, Cval, solinfo);
    end
end

%% Print and save table

colNames = strcat("deg_", string(degrees));

T = array2table(Bounds, 'VariableNames', cellstr(colNames));
T.Observable = string(obs(:,1));
T = movevars(T, 'Observable', 'Before', 1);

disp(T);

if saveTable
    writetable(T, 'lorenz_normalized_moment_bounds.csv');
    fprintf('Saved table: lorenz_normalized_moment_bounds.csv\n');
end

%% Defect plot for normalized z^4 bound

targetExp = [0 0 4];
targetD   = 12;
targetDegV = targetD - 1;

fprintf('\nComputing stored solution for z^4 defect plots...\n');

[Cz4, ~, defectPoly, vars] = solve_lorenz_moment_bound( ...
    targetExp, targetD, targetDegV, sigma, rho, beta, Kbox, opts);

fprintf('z^4 bound at relaxation degree %d: %.12g\n', targetD, Cz4);

if saveFigures
    plot_lorenz_defect(defectPoly, vars, sigma, rho, beta, Kbox, ...
        'xz', 'lorenz_z4_defect_xz.pdf');
    plot_lorenz_defect(defectPoly, vars, sigma, rho, beta, Kbox, ...
        'xy', 'lorenz_z4_defect_xy.pdf');
end

fprintf('\nFinished Lorenz time-average bound example.\n');

%% Local functions

function [Cval, solinfo, defect, vars] = solve_lorenz_moment_bound( ...
    expnt, d, degV, sigma, rho, beta, Kbox, opts)

    yalmip('clear');

    sdpvar X Y Z C
    vars = [X; Y; Z];

    x = Kbox.mid(1) + Kbox.rad(1)*X;
    y = Kbox.mid(2) + Kbox.rad(2)*Y;
    z = Kbox.mid(3) + Kbox.rad(3)*Z;

    fx = sigma*(y - x);
    fy = rho*x - y - x*z;
    fz = x*y - beta*z;

    fScaled = [
        fx/Kbox.rad(1);
        fy/Kbox.rad(2);
        fz/Kbox.rad(3)
    ];

    l = expnt(1);
    m = expnt(2);
    n = expnt(3);

    normFactor = beta^((l+m)/2) * (rho - 1)^((l+m)/2 + n);
    Phi = x^l * y^m * z^n / normFactor;

    [V, vcoeff] = symmetric_polynomial_lorenz(vars, degV);

    LfV = jacobian(V, vars)*fScaled;

    % Upper-bound defect:
    %     C - Phi - grad(V).f >= 0.
    defect = C - Phi - LfV;

    g = [
        1 - X^2;
        1 - Y^2;
        1 - Z^2
    ];

    constraints = [];

    [sigma0, c0] = sos_polynomial(vars, d, true);
    constraints = [constraints, sos(sigma0)];

    rhs = sigma0;
    allcoeff = [vcoeff; c0];

    for j = 1:length(g)

        [sigj, cj] = sos_polynomial(vars, d - 2, true);

        constraints = [constraints, sos(sigj)]; 
        rhs = rhs + sigj*g(j);
        allcoeff = [allcoeff; cj]; 
    end

    constraints = [constraints, coefficients(defect - rhs, vars) == 0];

    decvars = [C; allcoeff]; 

    sol = optimize(constraints, C, opts);

    Cval = value(C);
    solinfo = sol.info;

    if sol.problem ~= 0
        warning('Solver status for d=%d, observable [%d %d %d]: %s', ...
            d, expnt(1), expnt(2), expnt(3), sol.info);
    end
end

function [V, coeff] = symmetric_polynomial_lorenz(vars, degV)
% Polynomial invariant under (X,Y,Z) -> (-X,-Y,Z).

    X = vars(1);
    Y = vars(2);
    Z = vars(3);

    V = 0;
    coeff = [];

    for a = 0:degV
        for b = 0:(degV-a)
            for c = 0:(degV-a-b)

                if mod(a+b,2) == 0
                    q = sdpvar(1,1);
                    coeff = [coeff; q]; 
                    V = V + q*X^a*Y^b*Z^c;
                end
            end
        end
    end
end

function [p, coeffs] = sos_polynomial(vars, degMax, useSymmetry)

    if degMax < 0
        p = 0;
        coeffs = [];
        return
    end

    degMax = 2*floor(degMax/2);

    [p, coeffs] = polynomial_from_monomials(vars, degMax, useSymmetry);
end

function [p, coeffs] = polynomial_from_monomials(vars, degMax, useSymmetry)

    X = vars(1);
    Y = vars(2);
    Z = vars(3);

    p = 0;
    coeffs = [];

    for a = 0:degMax
        for b = 0:(degMax-a)
            for c = 0:(degMax-a-b)

                if useSymmetry && mod(a+b,2) ~= 0
                    continue
                end

                q = sdpvar(1,1);
                coeffs = [coeffs; q]; 
                p = p + q*X^a*Y^b*Z^c;
            end
        end
    end
end

function plot_lorenz_defect(defectPoly, vars, sigma, rho, beta, Kbox, plane, fileName)

    fOrig = @(~,u) [
        sigma*(u(2)-u(1));
        rho*u(1)-u(2)-u(1)*u(3);
        u(1)*u(2)-beta*u(3)
    ];

    optsODE = odeset('RelTol',1e-9,'AbsTol',1e-11);

    [~, u] = ode45(fOrig, [0 160], [1; 1; 1], optsODE);

    keep = round(0.35*size(u,1)):size(u,1);
    u = u(keep,:);

    Xa = (u(:,1) - Kbox.mid(1))/Kbox.rad(1);
    Ya = (u(:,2) - Kbox.mid(2))/Kbox.rad(2);
    Za = (u(:,3) - Kbox.mid(3))/Kbox.rad(3);

    Da = eval_yalmip_polynomial(defectPoly, vars, Xa, Ya, Za);
    Da = max(Da, 1e-14);
    logDa = log10(Da);

    fig = figure;
    set(fig,'Color','w','Units','centimeters','Position',[2 2 12 9]);

    hold on
    box on
    grid on

    switch lower(plane)

        case 'xz'
            scatter(u(:,1), u(:,3), 10, logDa, 'filled', ...
                'MarkerFaceAlpha', 1, ...
                'MarkerEdgeAlpha', 1);

            xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
            ylabel('$z$','Interpreter','latex','FontSize',24,'FontWeight','bold');

        case 'xy'
            scatter(u(:,1), u(:,2), 10, logDa, 'filled', ...
                'MarkerFaceAlpha', 1, ...
                'MarkerEdgeAlpha', 1);

            xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
            ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold');

        otherwise
            error('Unknown plane option: %s', plane);
    end

    colormap(parula);
    colorbar;

    set(gca,'FontSize',24,'TickLabelInterpreter','latex');

    save_pdf_figure(fig,fileName);
end

function val = eval_yalmip_polynomial(p, vars, X, Y, Z)

    [coef, monom] = coefficients(p, vars);
    coef = value(coef);

    val = zeros(size(X));

    for k = 1:length(coef)

        exps = zeros(1,3);

        for j = 1:3
            exps(j) = degree(monom(k), vars(j));
        end

        val = val + coef(k).*X.^exps(1).*Y.^exps(2).*Z.^exps(3);
    end

    val = real(val);
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