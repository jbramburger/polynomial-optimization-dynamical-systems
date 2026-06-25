%% rucklidge_pointwise_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Pointwise bounds on the Rucklidge attractor
%
% This script computes SOS pointwise upper and lower bounds on z for the
% Rucklidge system
%
%     x' = -a x + b y - y z,
%     y' = x,
%     z' = -z + y^2.
%
% The computation is performed in scaled coordinates on a ball containing a
% numerically observed attractor. The script sweeps over lambda values and
% SOS degrees, then plots projected upper-bound level sets in the (x,z)
% plane.
%
% Pressing Run should reproduce the numerical output and save:
%
%     rucklidge_lambda_sweep_z_bounds_ball.mat
%     rucklidge_bounds_lambda.pdf
%     rucklidge_upper_bound_levelsets_xz.pdf
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

a = 2;
b = 6.7;

lambda_vals = logspace(-2,2,101);
sos_degrees = [2 4 6 8];

lambda_plot = 1;
degrees_plot = [4 6 8];

solver_name = 'mosek';
verbose_solver = 0;

saveFigures = true;
saveData = true;

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.orange = [230, 120, 20]/255;
S.cyan   = [0, 170, 200]/255;
S.black  = [0, 0, 0];
S.gray   = [0.40, 0.40, 0.40];

plotColors = {S.gray,S.blue,S.red,S.green,S.orange,S.cyan};

%% Simulate attractor

f_ode = @(~,u) [
    -a*u(1) + b*u(2) - u(2)*u(3);
     u(1);
    -u(3) + u(2)^2
];

T = 1000;
u0 = [1; 0; 1];

opts_ode = odeset('RelTol',1e-9,'AbsTol',1e-11);

fprintf('\n============================================================\n');
fprintf('Rucklidge pointwise bounds\n');
fprintf('============================================================\n\n');

fprintf('Simulating Rucklidge attractor...\n');

[t,u] = ode45(f_ode,[0 T],u0,opts_ode); 

u_attr = u(round(0.5*length(t)):end,:);

zmin_sim = min(u_attr(:,3));
zmax_sim = max(u_attr(:,3));

%% Scaling and ball radius

safety = 1.10;

sx = safety*max(abs(u_attr(:,1)));
sy = safety*max(abs(u_attr(:,2)));
sz = safety*max(abs(u_attr(:,3)));

X_attr = u_attr(:,1)/sx;
Y_attr = u_attr(:,2)/sy;
Z_attr = u_attr(:,3)/sz;

Rball = 1.05*max(sqrt(X_attr.^2 + Y_attr.^2 + Z_attr.^2));

fprintf('Scaling values:\n');
fprintf('  sx    = %.8g\n',sx);
fprintf('  sy    = %.8g\n',sy);
fprintf('  sz    = %.8g\n',sz);
fprintf('  Rball = %.8g\n\n',Rball);

%% Scaled polynomial system

sdpvar X Y Z
vars = [X; Y; Z];

x = sx*X;
y = sy*Y;
z = sz*Z;

fX = (-a*x + b*y - y*z)/sx;
fY = x/sy;
fZ = (-z + y^2)/sz;

f = [fX; fY; fZ];

g = Rball^2 - X^2 - Y^2 - Z^2;

%% Solver options

solver_opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',1, ...
    'cachesolvers',1);

%% Part 1: lambda sweep

upper_bounds = nan(length(sos_degrees), length(lambda_vals));
lower_bounds = nan(length(sos_degrees), length(lambda_vals));

for id = 1:length(sos_degrees)

    d = sos_degrees(id);
    degV = d - 1;

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Lambda sweep: d = %d, deg(V) <= %d\n', d, degV);
    fprintf('------------------------------------------------------------\n');

    for ilam = 1:length(lambda_vals)

        lam = lambda_vals(ilam);

        fprintf('lambda = %.4e ... ', lam);

        [Cup, stat_up] = pointwise_bound_ball( ...
            vars, f, g, z, lam, d, degV, solver_opts);

        [Cminus, stat_low] = pointwise_bound_ball( ...
            vars, f, g, -z, lam, d, degV, solver_opts);

        if stat_up == 0
            upper_bounds(id, ilam) = Cup;
        end

        if stat_low == 0
            lower_bounds(id, ilam) = -Cminus;
        end

        fprintf('upper = %.8g, lower = %.8g\n', ...
            upper_bounds(id, ilam), lower_bounds(id, ilam));
    end
end

%% Remove trivial bounds from plotted curves

trivial_upper =  Rball*sz;
trivial_lower = -Rball*sz;

tol_trivial = 1e-6;

upper_plot = upper_bounds;
lower_plot = lower_bounds;

upper_plot(abs(upper_plot - trivial_upper) <= ...
    tol_trivial*max(1,abs(trivial_upper))) = NaN;

lower_plot(abs(lower_plot - trivial_lower) <= ...
    tol_trivial*max(1,abs(trivial_lower))) = NaN;

%% Plot lambda sweep

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 12 9]);

hold on
box on

yline(zmax_sim, ...
    'Color', S.black, ...
    'LineWidth', 2, ...
    'HandleVisibility','off');

yline(zmin_sim, ...
    'Color', S.black, ...
    'LineWidth', 2, ...
    'HandleVisibility','off');

for id = 1:length(sos_degrees)

    c = plotColors{id};

    semilogx(lambda_vals, upper_plot(id,:), '-', ...
        'Color', c, ...
        'LineWidth', 3.0, ...
        'DisplayName', sprintf('$d=%d$', sos_degrees(id)));

    semilogx(lambda_vals, lower_plot(id,:), '-', ...
        'Color', c, ...
        'LineWidth', 3.0, ...
        'HandleVisibility','off');
end

set(gca,'XScale','log');
xlim([min(lambda_vals), max(lambda_vals)]);

xlabel('$\lambda$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$z$ bounds','Interpreter','latex','FontSize',24,'FontWeight','bold');

legend('Location','west','Interpreter','latex');

set(gca,'FontSize',24,'TickLabelInterpreter','latex','Layer','top');
grid on
axis tight

if saveFigures
    save_pdf_figure(fig1,'rucklidge_bounds_lambda.pdf');
end

%% Save sweep results

results.lambda_vals = lambda_vals;
results.sos_degrees = sos_degrees;
results.upper_bounds = upper_bounds;
results.lower_bounds = lower_bounds;
results.upper_plot = upper_plot;
results.lower_plot = lower_plot;
results.simulated_z_min = zmin_sim;
results.simulated_z_max = zmax_sim;
results.trivial_upper = trivial_upper;
results.trivial_lower = trivial_lower;
results.scaling = struct('sx',sx,'sy',sy,'sz',sz, ...
    'safety',safety,'Rball',Rball);
results.parameters = struct('a',a,'b',b);

if saveData
    save('rucklidge_lambda_sweep_z_bounds_ball.mat','results');
    fprintf('Saved data: rucklidge_lambda_sweep_z_bounds_ball.mat\n');
end

%% Part 2: projected V = C boundaries for lambda = 1

nx = 450;
ny = 450;
nz = 450;

xGrid = linspace(-Rball*sx, Rball*sx, nx);
yGrid = linspace(-Rball*sy, Rball*sy, ny);
zGrid = linspace(-Rball*sz, Rball*sz, nz);

[Xg_plot,Zg_plot] = meshgrid(xGrid,zGrid);

upperProj = cell(length(degrees_plot),1);
upperC = nan(length(degrees_plot),1);

fprintf('\n------------------------------------------------------------\n');
fprintf('Projected V = C boundaries at lambda = %.3g\n', lambda_plot);
fprintf('------------------------------------------------------------\n');

for id = 1:length(degrees_plot)

    d = degrees_plot(id);
    degV = d - 1;

    fprintf('\nComputing upper-bound certificate for d = %d\n',d);

    [Cup, stat_up, Vup_data] = pointwise_certificate_ball( ...
        vars, f, g, z, lambda_plot, d, degV, solver_opts);

    if stat_up == 0

        upperC(id) = Cup;

        upperProj{id} = projected_level_function( ...
            Vup_data, Cup, xGrid, yGrid, zGrid, sx, sy, sz, Rball);

        fprintf('C = %.8g\n', Cup);

    else

        fprintf('Certificate failed for d = %d\n', d);
    end
end

%% Plot projected upper-bound trapping boundaries

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 12 9]);

hold on
box on

plot(u_attr(:,1), u_attr(:,3), '.', ...
    'Color', 0.70*S.gray + 0.30*[1 1 1], ...
    'MarkerSize', 2, ...
    'HandleVisibility','off');

for id = 1:length(degrees_plot)

    if ~isempty(upperProj{id})
        contour(Xg_plot,Zg_plot,upperProj{id},[0 0], ...
            'Color', plotColors{id+1}, ...
            'LineWidth', 3.0, ...
            'DisplayName', sprintf('$d=%d$',degrees_plot(id)));
    end
end

xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$z$','Interpreter','latex','FontSize',24,'FontWeight','bold');

legend('Location','best','Interpreter','latex');

set(gca,'FontSize',24,'TickLabelInterpreter','latex','Layer','top');
grid on

xlim([-12, 12]);
ylim([-0.1, 18]);

if saveFigures
    save_pdf_figure(fig2,'rucklidge_upper_bound_levelsets_xz.pdf');
end

fprintf('\nFinished Rucklidge pointwise-bound example.\n');

%% Local functions

function [Cval, status] = pointwise_bound_ball(vars, f, g, Phi, lambda, d, degV, solver_opts)

    sdpvar C

    [V, cV] = polynomial(vars, degV);

    [sig0, c_sig0] = polynomial(vars, d);
    [sig1, c_sig1] = polynomial(vars, d-2);

    [rho0, c_rho0] = polynomial(vars, d);
    [rho1, c_rho1] = polynomial(vars, d-2);

    trap_poly = C - V - lambda*jacobian(V,vars)*f;
    obs_poly  = V - Phi;

    constraints = [
        sos(sig0), ...
        sos(sig1), ...
        sos(rho0), ...
        sos(rho1), ...
        coefficients(trap_poly - sig0 - sig1*g, vars) == 0, ...
        coefficients(obs_poly  - rho0 - rho1*g, vars) == 0
        ];

    coeffs = [C; cV; c_sig0; c_sig1; c_rho0; c_rho1];

    sol = solvesos(constraints, C, solver_opts, coeffs);
    status = sol.problem;

    if status == 0
        Cval = value(C);
    else
        Cval = NaN;
    end
end

function [Cval, status, Vdata] = pointwise_certificate_ball( ...
    vars, f, g, Phi, lambda, d, degV, solver_opts)

    sdpvar C

    [V, cV] = polynomial(vars, degV);

    [sig0, c_sig0] = polynomial(vars, d);
    [sig1, c_sig1] = polynomial(vars, d-2);

    [rho0, c_rho0] = polynomial(vars, d);
    [rho1, c_rho1] = polynomial(vars, d-2);

    trap_poly = C - V - lambda*jacobian(V,vars)*f;
    obs_poly  = V - Phi;

    constraints = [
        sos(sig0), ...
        sos(sig1), ...
        sos(rho0), ...
        sos(rho1), ...
        coefficients(trap_poly - sig0 - sig1*g, vars) == 0, ...
        coefficients(obs_poly  - rho0 - rho1*g, vars) == 0
        ];

    coeffs = [C; cV; c_sig0; c_sig1; c_rho0; c_rho1];

    sol = solvesos(constraints, C, solver_opts, coeffs);
    status = sol.problem;

    if status == 0

        Cval = value(C);

        [Vcoeff, Vmonom] = coefficients(V,vars);

        Vdata.coeff = value(Vcoeff);
        Vdata.powers = monomial_powers(Vmonom, vars);

    else

        Cval = NaN;
        Vdata = [];
    end
end

function Fproj = projected_level_function( ...
    Vdata, Cval, xGrid, yGrid, zGrid, sx, sy, sz, Rball)

    nx = length(xGrid);
    nz = length(zGrid);

    Fproj = inf(nz,nx);

    Xs = xGrid/sx;
    Zs = zGrid/sz;

    [Xg,Zg] = meshgrid(Xs,Zs);

    for ky = 1:length(yGrid)

        Ys = yGrid(ky)/sy;

        inside = Xg.^2 + Ys.^2 + Zg.^2 <= Rball^2;

        Vval = eval_poly_scaled(Vdata.coeff, Vdata.powers, Xg, Ys, Zg);

        Ftemp = Vval - Cval;
        Ftemp(~inside) = inf;

        Fproj = min(Fproj, Ftemp);
    end

    outside_xz = Xg.^2 + Zg.^2 > Rball^2;
    Fproj(outside_xz) = NaN;
end

function val = eval_poly_scaled(coeff, powers, X, Y, Z)

    val = zeros(size(X));

    for k = 1:length(coeff)
        val = val + coeff(k) ...
            .* (X.^powers(k,1)) ...
            .* (Y.^powers(k,2)) ...
            .* (Z.^powers(k,3));
    end
end

function P = monomial_powers(monomials, vars)

    P = zeros(length(monomials), length(vars));

    for i = 1:length(monomials)
        for j = 1:length(vars)
            P(i,j) = degree(monomials(i), vars(j));
        end
    end
end

function save_pdf_figure(fig,outfilename)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',outfilename);
    fprintf('Saved figure: %s\n',outfilename);
end