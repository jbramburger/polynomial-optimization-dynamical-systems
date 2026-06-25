%% backward_vdp_roa_sos.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: Region-of-attraction estimates for the backward Van der Pol oscillator
%
% This script computes SOS inner approximations of the region of attraction
% for the origin in the backward-time Van der Pol oscillator
%
%     x' = -y,
%     y' = x - mu (1 - x^2) y.
%
% The certified sets have the form
%
%     { (x,y) : V(x,y) <= 1 },
%
% where V is a polynomial Lyapunov function. The quadratic part of V is
% optimized together with higher-order polynomial terms. For each degree,
% the script alternates between:
%
%   1. fixing V and finding an SOS multiplier, and
%   2. fixing the multiplier and improving V.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     backward_vdp_roa_c1_variableV.pdf
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

mu = 1;

degrees       = [4 6 8];
max_alt_iters = 8;

cval       = 1;
eps_sos    = 1e-7;
eps_update = 1e-9;

solver_name    = 'mosek';
verbose_solver = 0;
saveFigure     = true;

figureFile = 'backward_vdp_roa_c1_variableV.pdf';

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

region_colors = {S.blue, S.red, S.green, S.orange, S.purple, S.cyan};

%% Problem setup

sdpvar x y
z = [x; y];

f = [-y;
      x - mu*(1 - x^2)*y];

r2 = x^2 + y^2;
p1 = r2;
p2 = r2^2;

%% Plotting window and grid

xmin = -3; xmax = 3;
ymin = -4; ymax = 4;

xgrid = linspace(xmin,xmax,600);
ygrid = linspace(ymin,ymax,600);
[X,Y] = meshgrid(xgrid,ygrid);

dx = xgrid(2)-xgrid(1);
dy = ygrid(2)-ygrid(1);

%% Compute unstable limit cycle by integrating backward in time

odefun = @(~,z) [-z(2);
                  z(1) - mu*(1 - z(1)^2)*z(2)];

odeopts = odeset('RelTol',1e-10,'AbsTol',1e-12);

[~, zc] = ode45(odefun, [0 -120], [2;0], odeopts);

keep = round(0.55*size(zc,1)):size(zc,1);
zcycle = zc(keep,:);

% Samples slightly inside the observed basin boundary.
boundary_samples = 0.92*zcycle(1:10:end,:);

%% Vector field grid for plotting

nx = 15;
ny = 15;

[Xg,Yg] = meshgrid(linspace(xmin,xmax,nx), linspace(ymin,ymax,ny));

Ug = -Yg;
Vg = Xg - mu*(1 - Xg.^2).*Yg;

Spd = sqrt(Ug.^2 + Vg.^2);
Ug = Ug ./ max(Spd,1e-12);
Vg = Vg ./ max(Spd,1e-12);

%% Main SOS loop

fprintf('\n============================================================\n');
fprintf('Backward Van der Pol ROA inner approximations\n');
fprintf('============================================================\n\n');

fprintf('mu            = %.6g\n', mu);
fprintf('certified c   = %.6g\n', cval);
fprintf('degrees       = ');
fprintf('%d ', degrees);
fprintf('\n\n');

results = struct();

prev_P = [];
prev_exps = [];
prev_coeff = [];

for id = 1:length(degrees)

    d = degrees(id);

    fprintf('\n============================================================\n');
    fprintf('Degree d = %d\n', d);
    fprintf('============================================================\n');

    [exps, nfree] = even_exponents_2d_from_degree4(d);

    if isempty(prev_P)
        P_current = eye(2);
        coeff_current = zeros(nfree,1);
    else
        P_current = prev_P;
        coeff_current = warm_start_coefficients(exps, prev_exps, prev_coeff);
    end

    V_current = build_V_expr(P_current, exps, coeff_current, x, y);

    best_area = -inf;
    best_P = P_current;
    best_coeff = coeff_current;
    best_V_expr = V_current;

    for iter = 1:max_alt_iters

        fprintf('\n  Alternation %d\n', iter);

        %% Step 1: fixed V, find SOS multiplier

        V_fixed = clean(V_current, 1e-12);

        sdeg = d + 2;
        if mod(sdeg,2) == 1
            sdeg = sdeg + 1;
        end

        [feasible_s, s_current] = find_multiplier_fixed_V_c1( ...
            V_fixed, f, x, y, p1, p2, eps_sos, sdeg, opts);

        if ~feasible_s
            warning('Could not certify current V for degree %d.', d);
            break;
        end

        %% Estimate area of current certified region {V <= 1}

        Vvals = eval_V_grid(P_current, exps, coeff_current, X, Y);
        mask = Vvals <= cval;

        area_est = sum(mask(:))*dx*dy;

        fprintf('    Estimated certified area = %.6g\n', area_est);

        if area_est > best_area
            best_area = area_est;
            best_P = P_current;
            best_coeff = coeff_current;
            best_V_expr = V_current;
        end

        %% Step 2: fixed multiplier, improve V

        [V_candidate, P_candidate, coeff_candidate, success] = ...
            improve_V_fixed_c1_s( ...
                P_current, coeff_current, exps, f, x, y, p1, p2, ...
                eps_update, s_current, boundary_samples, opts);

        if success

            accepted = false;
            alphas = [1, 0.5, 0.25, 0.125, 0.0625];

            for ia = 1:length(alphas)

                alpha = alphas(ia);

                P_trial = (1-alpha)*P_current + alpha*P_candidate;
                coeff_trial = (1-alpha)*coeff_current + alpha*coeff_candidate;

                V_trial = build_V_expr(P_trial, exps, coeff_trial, x, y);
                V_trial = clean(V_trial, 1e-10);

                [feasible_trial, ~] = find_multiplier_fixed_V_c1( ...
                    V_trial, f, x, y, p1, p2, eps_sos, sdeg, opts);

                if feasible_trial
                    P_current = P_trial;
                    coeff_current = coeff_trial;
                    V_current = V_trial;
                    accepted = true;

                    fprintf('    Accepted V-update with alpha = %.4g\n', alpha);
                    break;
                end
            end

            if ~accepted
                fprintf('    Candidate V could not be re-certified; retaining previous V.\n');
            end

        else
            fprintf('    V-update failed; retaining previous V.\n');
        end
    end

    results(id).degree = d;
    results(id).c = cval;
    results(id).P = best_P;
    results(id).exps = exps;
    results(id).coeff = best_coeff;
    results(id).V_expr = best_V_expr;
    results(id).area = best_area;

    prev_P = best_P;
    prev_exps = exps;
    prev_coeff = best_coeff;
end

%% Plot certified regions

fig = figure;
set(fig,'Color','w','Units','centimeters','Position',[2 2 13 11]);
hold on; box on;

show_degrees = degrees;

quiver(Xg,Yg,Ug,Vg,0.30, ...
    'Color', 0.86*[1 1 1], ...
    'LineWidth', 0.6, ...
    'MaxHeadSize', 0.55);

h_roa = gobjects(length(show_degrees),1);
leg_entries = cell(length(show_degrees),1);

h_cycle = plot(zcycle(:,1), zcycle(:,2), '-', ...
    'Color', S.black, ...
    'LineWidth', 4);

for j = 1:length(show_degrees)

    dshow = show_degrees(j);
    rid = find([results.degree] == dshow, 1);

    if isempty(rid) || isnan(results(rid).area) || results(rid).area < 0
        continue;
    end

    Vvals = eval_V_grid(results(rid).P, results(rid).exps, ...
        results(rid).coeff, X, Y);

    contour(X,Y,Vvals,[1 1], ...
        'Color', region_colors{j}, ...
        'LineWidth', 3);

    h_roa(j) = plot(nan,nan,'-', ...
        'Color', region_colors{j}, ...
        'LineWidth', 3);

    leg_entries{j} = sprintf('$d=%d$', dshow);
end

plot(0,0,'o', ...
    'MarkerSize', 20, ...
    'MarkerFaceColor', S.green, ...
    'MarkerEdgeColor', S.black, ...
    'LineWidth', 1);

axis([-2.2 2.2 -3 3]);
axis equal
grid on

valid = isgraphics(h_roa);
h_leg = [h_roa(valid); h_cycle];
leg_labels = [leg_entries(valid); {'unstable limit cycle'}];

legend(h_leg, leg_labels, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'FontSize',20, ...
    'Box','on');

set(gca,'FontSize',24)
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold');

if saveFigure
    save_pdf_figure(fig,figureFile);
end

%% Print final results

fprintf('\n============================================================\n');
fprintf('Computed certified areas for fixed level c = 1\n');
fprintf('============================================================\n');

for id = 1:length(results)
    fprintf('  degree %2d: area = %.8g\n', ...
        results(id).degree, results(id).area);
end

fprintf('\nFinished backward Van der Pol ROA example.\n');

%% Local functions

function [exps, nfree] = even_exponents_2d_from_degree4(d)
% Free monomials of total degree 4,6,...,d.
% The quadratic part is represented separately by a variable matrix P.

    exps = [];

    for total = 4:2:d
        for i = 0:total
            j = total - i;
            exps = [exps; i j]; 
        end
    end

    nfree = size(exps,1);
end

function V = build_V_expr(P, exps, coeff, x, y)
% Build V(x,y) = z'Pz plus higher-order polynomial terms.

    z = [x; y];
    V = z.'*P*z;

    for k = 1:size(exps,1)
        V = V + coeff(k)*x^exps(k,1)*y^exps(k,2);
    end
end

function Vvals = eval_V_grid(P, exps, coeff, X, Y)
% Evaluate V on a grid.

    Vvals = P(1,1)*X.^2 + 2*P(1,2)*X.*Y + P(2,2)*Y.^2;

    for k = 1:size(exps,1)
        Vvals = Vvals + coeff(k)*X.^exps(k,1).*Y.^exps(k,2);
    end
end

function [isfeas, s_expr_value] = find_multiplier_fixed_V_c1( ...
    V_fixed, f, x, y, p1, p2, eps_cert, sdeg, opts)

    z = [x; y];

    [s, coeff_s] = even_polynomial_2d(x, y, sdeg);
    dVf = jacobian(V_fixed,z)*f;

    cons = [
        sos(s), ...
        sos(V_fixed - eps_cert*p1), ...
        sos(-dVf - eps_cert*p2 - s*(1 - V_fixed))
        ];

    sol = solvesos(cons, [], opts, coeff_s);

    isfeas = (sol.problem == 0);

    if isfeas
        s_expr_value = replace(s, coeff_s, value(coeff_s));
        s_expr_value = clean(s_expr_value, 1e-10);
    else
        s_expr_value = [];
    end
end

function [V_new_expr, P_val, coeff_val, success] = improve_V_fixed_c1_s( ...
    P_old, coeff_old, exps, f, x, y, p1, p2, eps_update, ...
    s_fixed, boundary_samples, opts)

    z = [x; y];

    nfree = size(exps,1);

    sdpvar p11 p12 p22
    P = [p11 p12; p12 p22];

    coeff = sdpvar(nfree,1);

    V = build_V_expr(P, exps, coeff, x, y);
    dVf = jacobian(V,z)*f;

    delta = 1e-6;
    P_step = 0.5;
    coeff_step = 10;

    cons = [];

    cons = [cons, P >= delta*eye(2)];
    cons = [cons, sos(V - eps_update*p1)];
    cons = [cons, sos(-dVf - eps_update*p2 - s_fixed*(1 - V))];

    cons = [cons, -P_step <= p11 - P_old(1,1) <= P_step];
    cons = [cons, -P_step <= p12 - P_old(1,2) <= P_step];
    cons = [cons, -P_step <= p22 - P_old(2,2) <= P_step];

    if nfree > 0
        cons = [cons, -coeff_step <= coeff - coeff_old <= coeff_step];
    end

    obj = 0;
    for k = 1:size(boundary_samples,1)
        obj = obj + replace(V, [x y], boundary_samples(k,:));
    end

    sol = solvesos(cons, obj, opts, [p11; p12; p22; coeff]);

    success = (sol.problem == 0);

    if ~success
        fprintf('    V-update solver status: %s\n', sol.info);
    end

    if success
        P_val = value(P);
        coeff_val = value(coeff);
    else
        P_val = P_old;
        coeff_val = coeff_old;
    end

    V_new_expr = build_V_expr(P_val, exps, coeff_val, x, y);
    V_new_expr = clean(V_new_expr, 1e-10);
end

function [p, coeff] = even_polynomial_2d(x, y, deg)
% Polynomial invariant under (x,y) -> (-x,-y).
% Keeps monomials x^i y^j with i+j even and i+j <= deg.

    coeff = [];
    p = 0;

    for total = 0:deg
        if mod(total,2) ~= 0
            continue;
        end

        for i = 0:total
            j = total - i;
            a = sdpvar(1,1);
            coeff = [coeff; a]; 
            p = p + a*x^i*y^j;
        end
    end
end

function coeff = warm_start_coefficients(exps, prev_exps, prev_coeff)
% Initialize coefficients using the previous degree solution.

    coeff = zeros(size(exps,1),1);

    if isempty(prev_exps) || isempty(prev_coeff)
        return;
    end

    for k = 1:size(exps,1)
        match = find(prev_exps(:,1) == exps(k,1) & ...
                     prev_exps(:,2) == exps(k,2), 1);

        if ~isempty(match)
            coeff(k) = prev_coeff(match);
        end
    end
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