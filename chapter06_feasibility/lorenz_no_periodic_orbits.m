%% lorenz_no_periodic_orbits.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: SOS exclusion of periodic orbits in the Lorenz equations
%
% This script certifies parameter ranges for which the Lorenz equations have
% no nontrivial periodic orbits inside a compact absorbing set.
%
% The certificate searches for a polynomial auxiliary function V satisfying
%
%     grad(V) * f <= - ||f||^2
%
% on a scaled absorbing ball K. This implies that every omega-limit set in K
% is contained in the equilibrium set, excluding nontrivial periodic orbits.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     lorenz_no_periodic_certification.pdf
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

sigma = 10;
beta  = 8/3;

rho_start = 0;
rho_step  = 0.01;
rho_max   = 14;

degrees = [1 3 5 7 9 11 13];

solver_name = 'mosek';
verbose_solver = 0;

use_conservative_radius = true;
R_factor_practical = 1.8;
extra_sdeg = 0;

saveFigure = true;
figureFile = 'lorenz_no_periodic_certification.pdf';

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',1, ...
    'cachesolvers',1);

%% Book plotting style

blue = [0, 92, 175]/255;

%% Symbolic scaled variables

sdpvar X Y Z
gK = 1 - (X^2 + Y^2 + Z^2);

rho_values = rho_start:rho_step:rho_max;

certified_degree = NaN(size(rho_values));
certified_status = strings(size(rho_values));
solver_info = strings(size(rho_values));

current_degree_index = 1;

fprintf('\n============================================================\n');
fprintf('SOS exclusion of periodic orbits for the Lorenz equations\n');
fprintf('============================================================\n\n');

fprintf('sigma   = %.8g\n', sigma);
fprintf('beta    = %.8g\n', beta);
fprintf('rho grid: %.4g : %.4g : %.4g\n', rho_start, rho_step, rho_max);
fprintf('degrees : ');
fprintf('%d ', degrees);
fprintf('\n\n');

%% Main continuation loop in rho

for irho = 1:length(rho_values)

    rho = rho_values(irho);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('rho = %.4f\n', rho);
    fprintf('------------------------------------------------------------\n');

    certified_here = false;

    while current_degree_index <= length(degrees)

        dV = degrees(current_degree_index);

        fprintf('  Trying degree d = %d ... ', dV);

        [isfeas, solinfo] = certify_lorenz_no_periodic( ...
            sigma, beta, rho, dV, X, Y, Z, gK, ...
            use_conservative_radius, R_factor_practical, extra_sdeg, opts);

        if isfeas

            fprintf('certified.\n');

            certified_here = true;
            certified_degree(irho) = dV;
            certified_status(irho) = "certified";
            solver_info(irho) = string(solinfo);

            break;

        else

            fprintf('failed (%s).\n', solinfo);

            current_degree_index = current_degree_index + 1;

            if current_degree_index <= length(degrees)
                fprintf('  Increasing degree to d = %d.\n', ...
                    degrees(current_degree_index));
            end
        end
    end

    if ~certified_here
        fprintf('\nNo certification found up to maximum degree d = %d.\n', degrees(end));
        fprintf('Stopping continuation at rho = %.4f.\n', rho);

        certified_status(irho) = "failed";
        solver_info(irho) = "failed all degrees";

        break;
    end
end

%% Trim results

last_valid = find(~isnan(certified_degree), 1, 'last');

if isempty(last_valid)
    warning('No rho values certified.');
    last_valid = 0;
end

rho_cert = rho_values(1:last_valid);
deg_cert = certified_degree(1:last_valid);

%% Print certified ranges

fprintf('\n============================================================\n');
fprintf('Certified ranges by degree\n');
fprintf('============================================================\n');

if last_valid == 0

    fprintf('No certified rho values found.\n');

else

    for jd = 1:length(degrees)

        d = degrees(jd);
        idx = find(deg_cert == d);

        if isempty(idx)
            continue;
        end

        intervals = contiguous_intervals(idx);

        for q = 1:size(intervals,1)

            i1 = intervals(q,1);
            i2 = intervals(q,2);

            fprintf('  degree %2d: rho in [%.4f, %.4f]\n', ...
                d, rho_cert(i1), rho_cert(i2));
        end
    end
end

if last_valid < length(rho_values)
    fprintf('\nCertification stopped before rho = %.4f.\n', rho_max);
    if last_valid > 0
        fprintf('Last certified rho = %.4f.\n', rho_values(last_valid));
    end
end

%% Store numerical results

results = struct();
results.sigma = sigma;
results.beta = beta;
results.rho_values = rho_values;
results.certified_degree = certified_degree;
results.certified_status = certified_status;
results.solver_info = solver_info;
results.degrees = degrees;
results.rho_step = rho_step;
results.rho_max = rho_max;

%% Plot horizontal certification bars

if last_valid > 0

    fig = figure;
    set(fig,'Color','w','Units','centimeters','Position',[2 2 13 7]);

    hold on
    box on

    bar_height = 0.65;

    for jd = 1:length(degrees)

        d = degrees(jd);
        idx = find(deg_cert == d);

        if isempty(idx)
            continue;
        end

        intervals = contiguous_intervals(idx);

        for q = 1:size(intervals,1)

            i1 = intervals(q,1);
            i2 = intervals(q,2);

            r1 = rho_cert(i1);
            r2 = rho_cert(i2);

            patch([r1 r2 r2 r1], ...
                  [d-bar_height/2 d-bar_height/2 ...
                   d+bar_height/2 d+bar_height/2], ...
                  blue, ...
                  'FaceAlpha',1, ...
                  'EdgeColor','k', ...
                  'LineWidth',2);
        end
    end

    xlabel('$\rho$','Interpreter','latex','FontSize',18);
    ylabel('degree of $V$','Interpreter','latex','FontSize',18);

    yticks(degrees);
    ylim([min(degrees)-1 max(degrees)+1]);
    xlim([rho_start rho_cert(end)]);

    set(gca, ...
        'FontSize',16, ...
        'LineWidth',0.8, ...
        'TickLabelInterpreter','latex', ...
        'Layer','top');

    grid on

    if saveFigure
        save_pdf_figure(fig,figureFile);
    end
end

fprintf('\nFinished Lorenz no-periodic-orbit certification example.\n');

%% Local functions

function [isfeas, solinfo] = certify_lorenz_no_periodic( ...
    sigma, beta, rho, dV, X, Y, Z, gK, ...
    use_conservative_radius, R_factor_practical, extra_sdeg, opts)

    vars = [X; Y; Z];

    %% Scaled absorbing set

    zc = rho + sigma;

    if use_conservative_radius

        safety = 1.05;
        R = safety * beta/(2*sqrt(beta - 1)) * (rho + sigma);

    else

        R = R_factor_practical*(rho + sigma);
    end

    Lx = R;
    Ly = R;
    Lz = R;

    %% Physical variables in terms of scaled variables

    x = Lx*X;
    y = Ly*Y;
    z = zc + Lz*Z;

    %% Lorenz vector field in physical variables

    f_phys = [
        sigma*(y - x);
        rho*x - y - x*z;
        x*y - beta*z
    ];

    %% Scaled vector field

    F = [
        f_phys(1)/Lx;
        f_phys(2)/Ly;
        f_phys(3)/Lz
    ];

    Fnorm2 = F.'*F;

    %% Symmetry-reduced auxiliary polynomial V

    [V, coeffV] = symmetric_polynomial_lorenz(X, Y, Z, dV);

    dVf = jacobian(V,vars)*F;

    q = -dVf - Fnorm2;

    %% Multiplier degree

    qdeg = max(dV + 1, 4);

    if mod(qdeg,2) == 1
        qdeg = qdeg + 1;
    end

    sdeg = qdeg - 2 + extra_sdeg;

    if sdeg < 0
        sdeg = 0;
    end

    if mod(sdeg,2) == 1
        sdeg = sdeg + 1;
    end

    %% SOS multiplier polynomial

    [sK, coeffS] = full_polynomial_3d(X, Y, Z, sdeg);

    %% SOS constraints

    cons = [
        sos(sK), ...
        sos(q - sK*gK)
        ];

    decision_vars = [coeffV; coeffS];

    try
        sol = solvesos(cons, [], opts, decision_vars);
        isfeas = (sol.problem == 0);
        solinfo = sol.info;
    catch ME
        isfeas = false;
        solinfo = ME.message;
    end
end

function [V, coeff] = symmetric_polynomial_lorenz(X, Y, Z, deg)
% Polynomial invariant under (X,Y,Z) -> (-X,-Y,Z).
% Keeps monomials X^i Y^j Z^k with i+j even and total degree <= deg.
% The constant term is omitted because it has no effect on grad(V)*f.

    coeff = [];
    V = 0;

    for total = 1:deg
        for i = 0:total
            for j = 0:(total-i)

                k = total - i - j;

                if mod(i+j,2) ~= 0
                    continue;
                end

                a = sdpvar(1,1);
                coeff = [coeff; a]; %#ok<AGROW>
                V = V + a*X^i*Y^j*Z^k;
            end
        end
    end
end

function [p, coeff] = full_polynomial_3d(X, Y, Z, deg)
% Full polynomial in X,Y,Z of total degree <= deg.
% This avoids occasional YALMIP orientation issues with polynomial(vars,deg).

    coeff = [];
    p = 0;

    for total = 0:deg
        for i = 0:total
            for j = 0:(total-i)

                k = total - i - j;

                a = sdpvar(1,1);
                coeff = [coeff; a]; %#ok<AGROW>
                p = p + a*X^i*Y^j*Z^k;
            end
        end
    end
end

function intervals = contiguous_intervals(idx)
% Convert sorted index list into contiguous intervals.

    if isempty(idx)
        intervals = [];
        return;
    end

    idx = idx(:);
    breaks = [1; find(diff(idx) > 1) + 1; length(idx)+1];

    intervals = zeros(length(breaks)-1,2);

    for k = 1:length(breaks)-1
        intervals(k,1) = idx(breaks(k));
        intervals(k,2) = idx(breaks(k+1)-1);
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