%% perturbed_hopf_3d_roa_sos.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: Region-of-attraction estimates for a perturbed 3D Hopf-type system
%
% This script computes SOS inner approximations of the region of attraction
% for the origin in the perturbed Hopf-type system
%
%     x' = -y + x(1-r^2)(r^2-rho^2) + eps*x*z,
%     y' =  x + y(1-r^2)(r^2-rho^2) + eps*y*z,
%     z' = -z + eta*x^2,
%
% where r^2 = x^2 + y^2.
%
% The certified sets have the form
%
%     { (x,y,z) : V(x,y,z) <= 1 },
%
% where V is a polynomial Lyapunov function. The basin boundary is first
% estimated numerically by ray bisection. These samples are not part of the
% certificate; they are used only to guide the objective when improving V.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     perturbed_hopf_3d_roa_c1.png
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

rho      = 0.60;
eps_coup = 0.15;
eta      = 0.25;

degrees       = [2 4 6];
max_alt_iters = 12;

cval       = 1;
eps_sos    = 1e-9;
eps_update = 1e-10;

solver_name    = 'mosek';
verbose_solver = 0;

saveFigure = true;
figureFile = 'perturbed_hopf_3d_roa_c1.png';

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',1);

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.purple = [170, 90, 160]/255;
S.cyan   = [0, 170, 200]/255;
S.black  = [0, 0, 0];

region_colors = {S.blue, S.red, S.green, S.purple, S.cyan};

%% Symbolic variables and vector field

sdpvar x y z
Xsym = [x; y; z];

r2 = x^2 + y^2;

f = [
    -y + x*(1-r2)*(r2-rho^2) + eps_coup*x*z;
     x + y*(1-r2)*(r2-rho^2) + eps_coup*y*z;
    -z + eta*x^2
];

r3sq = x^2 + y^2 + z^2;

p1 = r3sq;
p2 = r3sq^2;

%% Numerical vector field

odefun = @(~,u) [
    -u(2) + u(1)*(1-(u(1)^2+u(2)^2))*((u(1)^2+u(2)^2)-rho^2) ...
        + eps_coup*u(1)*u(3);
     u(1) + u(2)*(1-(u(1)^2+u(2)^2))*((u(1)^2+u(2)^2)-rho^2) ...
        + eps_coup*u(2)*u(3);
    -u(3) + eta*u(1)^2
];

odeopts = odeset('RelTol',1e-9,'AbsTol',1e-11);

%% Estimate basin boundary by ray bisection

fprintf('\n============================================================\n');
fprintf('Perturbed Hopf 3D ROA example\n');
fprintf('============================================================\n\n');

fprintf('Parameters:\n');
fprintf('  rho      = %.6g\n', rho);
fprintf('  eps_coup = %.6g\n', eps_coup);
fprintf('  eta      = %.6g\n\n', eta);

fprintf('Estimating basin-boundary samples...\n');

num_dirs     = 400;
directions   = fibonacci_sphere(num_dirs);

R_low        = 0.05;
R_high       = 3.0;
Tclass       = 120;
conv_tol     = 0.08;
bisect_steps = 25;

boundary_samples = zeros(num_dirs,3);

for k = 1:num_dirs

    dir = directions(k,:).';

    lo = R_low;
    hi = R_high;

    for j = 1:bisect_steps

        mid = 0.5*(lo+hi);
        u0 = mid*dir;

        if converges_to_origin(odefun, u0, Tclass, conv_tol, odeopts)
            lo = mid;
        else
            hi = mid;
        end
    end

    boundary_samples(k,:) = (lo*dir).';
end

alpha_inside = 0.97;
target_samples = alpha_inside*boundary_samples;

fprintf('Finished estimating %d boundary samples.\n\n', num_dirs);

%% Plotting grid for isosurfaces and volume estimates

xmin = -1.4; xmax = 1.4;
ymin = -1.4; ymax = 1.4;
zmin = -0.4; zmax = 1.0;

ngrid = 55;

xgrid = linspace(xmin,xmax,ngrid);
ygrid = linspace(ymin,ymax,ngrid);
zgrid = linspace(zmin,zmax,ngrid);

[X,Y,Z] = meshgrid(xgrid,ygrid,zgrid);

dvol = (xgrid(2)-xgrid(1)) ...
     * (ygrid(2)-ygrid(1)) ...
     * (zgrid(2)-zgrid(1));

%% Main SOS loop

results = struct();

prev_P = [];
prev_exps = [];
prev_coeff = [];

fprintf('Starting alternating SOS search.\n');

for id = 1:length(degrees)

    d = degrees(id);

    fprintf('\n============================================================\n');
    fprintf('Degree d = %d\n', d);
    fprintf('============================================================\n');

    [exps, nfree] = monomial_exponents_3d(3,d);

    if isempty(prev_P)
        R_init = 0.30;
        P_current = (1/R_init^2)*eye(3);
        coeff_current = zeros(nfree,1);
    else
        P_current = prev_P;
        coeff_current = warm_start_coefficients_3d(exps, prev_exps, prev_coeff);
    end

    V_current = build_V_expr_3d(P_current, exps, coeff_current, x, y, z);

    best_volume = -inf;
    best_V_expr = V_current;
    best_P = P_current;
    best_coeff = coeff_current;

    for iter = 1:max_alt_iters

        fprintf('\n  Alternation %d\n', iter);

        %% Step 1: fixed V, find an SOS multiplier

        V_fixed = clean(V_current,1e-12);

        sdeg = d + 2;
        if mod(sdeg,2) == 1
            sdeg = sdeg + 1;
        end

        [feasible_s, s_current] = find_multiplier_fixed_V_c1_3d( ...
            V_fixed, f, x, y, z, p1, p2, eps_sos, sdeg, opts);

        if ~feasible_s
            warning('Could not certify current V for degree %d.', d);
            break;
        end

        %% Estimate volume of certified region {V <= 1}

        Vvals = eval_V_grid_3d(P_current, exps, coeff_current, X, Y, Z);
        mask = Vvals <= cval;
        volume_est = sum(mask(:))*dvol;

        fprintf('    Estimated certified volume = %.6g\n', volume_est);

        if volume_est > best_volume
            best_volume = volume_est;
            best_V_expr = V_current;
            best_P = P_current;
            best_coeff = coeff_current;
        end

        %% Step 2: fixed multiplier, improve V

        if iter < max_alt_iters

            [V_candidate, P_candidate, coeff_candidate, success] = ...
                improve_V_fixed_c1_s_3d( ...
                    P_current, coeff_current, exps, f, x, y, z, ...
                    p1, p2, eps_update, s_current, target_samples, opts);

            if success

                accepted = false;
                theta_values = [1, 0.5, 0.25, 0.125, 0.0625];

                for itheta = 1:length(theta_values)

                    theta = theta_values(itheta);

                    P_trial = (1-theta)*P_current + theta*P_candidate;
                    coeff_trial = (1-theta)*coeff_current + theta*coeff_candidate;

                    V_trial = build_V_expr_3d(P_trial, exps, coeff_trial, x, y, z);
                    V_trial = clean(V_trial,1e-10);

                    [feasible_trial, ~] = find_multiplier_fixed_V_c1_3d( ...
                        V_trial, f, x, y, z, p1, p2, eps_update, sdeg, opts);

                    if feasible_trial
                        P_current = P_trial;
                        coeff_current = coeff_trial;
                        V_current = V_trial;

                        accepted = true;
                        fprintf('    Accepted V-update with theta = %.4g\n', theta);
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
    end

    results(id).degree = d;
    results(id).c = cval;
    results(id).P = best_P;
    results(id).exps = exps;
    results(id).coeff = best_coeff;
    results(id).V_expr = best_V_expr;
    results(id).volume = best_volume;

    prev_P = best_P;
    prev_exps = exps;
    prev_coeff = best_coeff;
end

%% Plot 3D certified regions and boundary samples

fig = figure;
set(fig,'Color','w','Units','centimeters','Position',[2 2 13 11]);

hold on
box on

% Certified sublevel surfaces, plotted from largest degree to smallest.
for id = length(results):-1:1

    if isnan(results(id).volume) || results(id).volume < 0
        continue;
    end

    Vvals = eval_V_grid_3d(results(id).P, results(id).exps, ...
        results(id).coeff, X, Y, Z);

    fv = isosurface(X,Y,Z,Vvals,cval);

    if isempty(fv.vertices)
        continue;
    end

    patch(fv, ...
        'FaceColor', region_colors{id}, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.10);
end

% Boundary surfaces as edges.
for id = 1:length(results)

    if isnan(results(id).volume) || results(id).volume < 0
        continue;
    end

    Vvals = eval_V_grid_3d(results(id).P, results(id).exps, ...
        results(id).coeff, X, Y, Z);

    fv = isosurface(X,Y,Z,Vvals,cval);

    if isempty(fv.vertices)
        continue;
    end

    patch(fv, ...
        'FaceColor', 'none', ...
        'EdgeColor', region_colors{id}, ...
        'LineWidth', 0.4, ...
        'EdgeAlpha', 0.3);
end

plot3(0,0,0,'o', ...
    'MarkerSize',9, ...
    'MarkerFaceColor',S.green, ...
    'MarkerEdgeColor',S.black, ...
    'LineWidth',0.9);

scatter3(boundary_samples(:,1), boundary_samples(:,2), boundary_samples(:,3), ...
    30, S.black, 'filled', ...
    'MarkerFaceAlpha',1, ...
    'MarkerEdgeAlpha',1);

axis([xmin xmax ymin ymax zmin zmax]);
axis equal
view(42,22);

xlabel('$x$','Interpreter','latex','FontSize',16);
ylabel('$y$','Interpreter','latex','FontSize',16);
zlabel('$z$','Interpreter','latex','FontSize',16);

set(gca, ...
    'FontSize',16, ...
    'LineWidth',2, ...
    'TickLabelInterpreter','latex', ...
    'Layer','top');

camlight headlight
grid on

h = gobjects(length(results)+1,1);
labels = cell(length(results)+1,1);

for id = 1:length(results)
    h(id) = patch(nan,nan,nan,region_colors{id}, ...
        'FaceAlpha',0.75, ...
        'EdgeColor','none');
    labels{id} = sprintf('$d=%d$', results(id).degree);
end

h(end) = scatter3(nan,nan,nan,30,S.black,'filled');
labels{end} = 'boundary samples';

legend(h, labels, ...
    'Interpreter','latex', ...
    'Location','northeast', ...
    'FontSize',16, ...
    'Box','on');

set(fig,'Renderer','opengl');

if saveFigure
    exportgraphics(fig, figureFile, ...
        'Resolution',600, ...
        'BackgroundColor','white');
    fprintf('Saved figure: %s\n', figureFile);
end

%% Print final results

fprintf('\n============================================================\n');
fprintf('Computed certified volumes for fixed level c = 1\n');
fprintf('============================================================\n');

for id = 1:length(results)
    fprintf('  degree %2d: volume = %.8g\n', ...
        results(id).degree, results(id).volume);
end

fprintf('\nFinished perturbed Hopf 3D ROA example.\n');

%% Local functions

function tf = converges_to_origin(odefun, u0, T, tol, odeopts)

    [~,u] = ode45(odefun,[0 T],u0,odeopts);
    tf = norm(u(end,:)) < tol;
end

function dirs = fibonacci_sphere(N)
% Approximately uniform directions on the sphere.

    dirs = zeros(N,3);
    phi = pi*(3 - sqrt(5));

    for k = 1:N
        y = 1 - 2*(k-1)/(N-1);
        r = sqrt(max(0,1-y^2));
        theta = phi*(k-1);

        dirs(k,:) = [cos(theta)*r, sin(theta)*r, y];
    end
end

function [exps, nfree] = monomial_exponents_3d(min_degree, max_degree)
% All monomials x^i y^j z^k with min_degree <= i+j+k <= max_degree.

    exps = [];

    if max_degree < min_degree
        nfree = 0;
        return;
    end

    for total = min_degree:max_degree
        for i = 0:total
            for j = 0:(total-i)
                k = total - i - j;
                exps = [exps; i j k]; 
            end
        end
    end

    nfree = size(exps,1);
end

function V = build_V_expr_3d(P, exps, coeff, x, y, z)
% Build V(x,y,z) = X'PX plus free higher-order polynomial terms.

    X = [x; y; z];
    V = X.'*P*X;

    for k = 1:size(exps,1)
        V = V + coeff(k)*x^exps(k,1)*y^exps(k,2)*z^exps(k,3);
    end
end

function Vvals = eval_V_grid_3d(P, exps, coeff, X, Y, Z)
% Evaluate V on a 3D grid.

    Vvals = P(1,1)*X.^2 + P(2,2)*Y.^2 + P(3,3)*Z.^2 ...
        + 2*P(1,2)*X.*Y + 2*P(1,3)*X.*Z + 2*P(2,3)*Y.*Z;

    for k = 1:size(exps,1)
        Vvals = Vvals ...
            + coeff(k)*X.^exps(k,1).*Y.^exps(k,2).*Z.^exps(k,3);
    end
end

function [isfeas, s_expr_value] = find_multiplier_fixed_V_c1_3d( ...
    V_fixed, f, x, y, z, p1, p2, eps_cert, sdeg, opts)

    X = [x; y; z];

    [s, coeff_s] = polynomial(X, sdeg);
    dVf = jacobian(V_fixed,X)*f;

    cons = [
        sos(s), ...
        sos(V_fixed - eps_cert*p1), ...
        sos(-dVf - eps_cert*p2 - s*(1 - V_fixed))
        ];

    sol = solvesos(cons, [], opts, coeff_s);

    isfeas = (sol.problem == 0);

    if isfeas
        s_expr_value = replace(s, coeff_s, value(coeff_s));
        s_expr_value = clean(s_expr_value,1e-10);
    else
        s_expr_value = [];
    end
end

function [V_new_expr, P_val, coeff_val, success] = improve_V_fixed_c1_s_3d( ...
    P_old, coeff_old, exps, f, x, y, z, p1, p2, eps_update, ...
    s_fixed, target_samples, opts)

    X = [x; y; z];
    nfree = size(exps,1);

    sdpvar p11 p12 p13 p22 p23 p33

    P = [p11 p12 p13;
         p12 p22 p23;
         p13 p23 p33];

    coeff = sdpvar(nfree,1);

    V = build_V_expr_3d(P, exps, coeff, x, y, z);
    dVf = jacobian(V,X)*f;

    delta = 1e-6;
    P_step = 0.5;
    coeff_step = 10;

    cons = [];

    cons = [cons, P >= delta*eye(3)];
    cons = [cons, sos(V - eps_update*p1)];
    cons = [cons, sos(-dVf - eps_update*p2 - s_fixed*(1 - V))];

    cons = [cons, -P_step <= p11 - P_old(1,1) <= P_step];
    cons = [cons, -P_step <= p12 - P_old(1,2) <= P_step];
    cons = [cons, -P_step <= p13 - P_old(1,3) <= P_step];
    cons = [cons, -P_step <= p22 - P_old(2,2) <= P_step];
    cons = [cons, -P_step <= p23 - P_old(2,3) <= P_step];
    cons = [cons, -P_step <= p33 - P_old(3,3) <= P_step];

    if nfree > 0
        cons = [cons, -coeff_step <= coeff - coeff_old <= coeff_step];
    end

    obj = 0;
    for k = 1:size(target_samples,1)
        obj = obj + replace(V, [x y z], target_samples(k,:));
    end

    decision_vars = [p11; p12; p13; p22; p23; p33; coeff];

    sol = solvesos(cons, obj, opts, decision_vars);

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

    V_new_expr = build_V_expr_3d(P_val, exps, coeff_val, x, y, z);
    V_new_expr = clean(V_new_expr,1e-10);
end

function coeff = warm_start_coefficients_3d(exps, prev_exps, prev_coeff)
% Initialize coefficients for the current degree using the previous degree.

    coeff = zeros(size(exps,1),1);

    if isempty(prev_exps) || isempty(prev_coeff)
        return;
    end

    for k = 1:size(exps,1)

        match = find(prev_exps(:,1) == exps(k,1) & ...
                     prev_exps(:,2) == exps(k,2) & ...
                     prev_exps(:,3) == exps(k,3), 1);

        if ~isempty(match)
            coeff(k) = prev_coeff(match);
        end
    end
end