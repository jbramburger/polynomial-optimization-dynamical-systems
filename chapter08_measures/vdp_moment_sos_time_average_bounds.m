%% vdp_moment_sos_time_average_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Invariant Measures, Ergodic Optimization, and Duality
% Example: Moment-SDP and SOS bounds for Van der Pol time averages
%
% This script compares moment-SDP bounds and SOS auxiliary-function bounds
% for the long-time average of
%
%     Phi(x,y) = x^2 + y^2
%
% over invariant probability measures of the Van der Pol oscillator
%
%     x' = y,
%     y' = (1 - x^2)y - x.
%
% The computation is performed in scaled variables
%
%     x = R X,   y = R Y,
%
% on the unit disk X^2 + Y^2 <= 1. The observable used in the SDP is
%
%     Phi(X,Y) = X^2 + Y^2,
%
% and the final reported bounds are multiplied by R^2.
%
% Pressing Run should reproduce the numerical output and save:
%
%     vdp_time_average_bounds.csv
%     vdp_time_average_defect.pdf
%     vdp_optimal_measure.pdf
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
R  = 4;

degrees = 2:2:20;
plotDegree = 20;

gridSize = 512;

solver_name = 'mosek';
verbose_solver = 0;

saveFigures = true;
saveTable = true;

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.black  = [0, 0, 0];

%% Solver options

optsMoment = sdpsettings('solver',solver_name,'verbose',verbose_solver);
optsSOS    = sdpsettings('solver',solver_name,'verbose',verbose_solver);

%% Storage

nDeg = numel(degrees);

momentLower = nan(nDeg,1);
momentUpper = nan(nDeg,1);
sosLower    = nan(nDeg,1);
sosUpper    = nan(nDeg,1);

upperMomentStore = cell(nDeg,1);
momentBasisStore = cell(nDeg,1);

upperDefectStore = cell(nDeg,1);
defectVarsStore  = cell(nDeg,1);

fprintf('\n============================================================\n');
fprintf('Van der Pol time-average bounds via moments and SOS\n');
fprintf('============================================================\n\n');

fprintf('mu              = %.6g\n', mu);
fprintf('scaling radius  = %.6g\n', R);
fprintf('degrees         = ');
fprintf('%d ', degrees);
fprintf('\n\n');

%% Main hierarchy loop

for kk = 1:nDeg

    d = degrees(kk);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Relaxation degree d = %d\n', d);
    fprintf('------------------------------------------------------------\n');

    %% Moment SDP bounds

    try

        [momentUpper(kk), yUpper, basisY] = solve_moment_sdp_vdp( ...
            d, R, optsMoment, +1);

        [momentLower(kk), ~, ~] = solve_moment_sdp_vdp( ...
            d, R, optsMoment, -1);

        upperMomentStore{kk} = yUpper;
        momentBasisStore{kk} = basisY;

        fprintf('Moment lower/upper: %.10f   %.10f\n', ...
            R^2*momentLower(kk), R^2*momentUpper(kk));

    catch ME

        warning('Moment SDP failed at degree %d:\n%s', d, ME.message);
    end

    %% SOS auxiliary-function bounds

    try

        [sosUpper(kk), defectUpper, defectVars] = solve_sos_bound_vdp( ...
            d, R, optsSOS, +1);

        [sosLower(kk), ~, ~] = solve_sos_bound_vdp( ...
            d, R, optsSOS, -1);

        upperDefectStore{kk} = defectUpper;
        defectVarsStore{kk}  = defectVars;

        fprintf('SOS lower/upper:    %.10f   %.10f\n', ...
            R^2*sosLower(kk), R^2*sosUpper(kk));

    catch ME

        warning('SOS problem failed at degree %d:\n%s', d, ME.message);
    end
end

%% Summary table

Results = table( ...
    degrees(:), ...
    R^2*momentLower, ...
    R^2*momentUpper, ...
    R^2*sosLower, ...
    R^2*sosUpper, ...
    'VariableNames', {'Degree','MomentLower','MomentUpper','SOSLower','SOSUpper'});

fprintf('\n============================================================\n');
fprintf('Summary of bounds in original variables\n');
fprintf('============================================================\n\n');

disp(Results);

if saveTable
    writetable(Results,'vdp_time_average_bounds.csv');
    fprintf('Saved table: vdp_time_average_bounds.csv\n');
end

%% Reconstruct optimal measure and upper-bound defect

kkPlot = find(degrees == plotDegree, 1);

if isempty(kkPlot)
    error('Requested plotDegree = %d is not included in degrees.', plotDegree);
end

if isempty(upperMomentStore{kkPlot}) || isempty(upperDefectStore{kkPlot})
    error('Degree %d did not solve successfully, so plots cannot be generated.', plotDegree);
end

fprintf('\nReconstructing optimal measure and defect at degree %d...\n', plotDegree);

yUpper = upperMomentStore{kkPlot};
basisY = momentBasisStore{kkPlot};

[densityCoeff, densityBasis] = reconstruct_density_from_moments( ...
    yUpper, basisY, plotDegree);

defectPolynomial = upperDefectStore{kkPlot};
defectVars = defectVarsStore{kkPlot};

%% Evaluate on grid

xx = linspace(-1,1,gridSize);
yy = linspace(-1,1,gridSize);

[Xg,Yg] = meshgrid(xx,yy);
mask = Xg.^2 + Yg.^2 <= 1;

densityGrid = nan(size(Xg));
defectGrid  = nan(size(Xg));

for q = 1:numel(Xg)

    if mask(q)

        densityGrid(q) = eval_poly_from_coeffs( ...
            densityCoeff, densityBasis, Xg(q), Yg(q));

        defectGrid(q) = value(replace( ...
            defectPolynomial, defectVars, [Xg(q); Yg(q)]));
    end
end

defectGrid = max(defectGrid,1e-14);

%% Figure 1: SOS defect

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 11 9]);

contourf(R*Xg,R*Yg,log10(defectGrid),35,'LineColor','none');
colormap(parula);
colorbar;

set(gca,'FontSize',42,'TickLabelInterpreter','latex');

xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

axis equal tight
box on

if saveFigures
    save_pdf_figure(fig1,'vdp_time_average_defect.pdf');
end

%% Figure 2: reconstructed optimal measure density

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 11 9]);

contourf(R*Xg,R*Yg,densityGrid,35,'LineColor','none');
colormap(parula);
colorbar;

set(gca,'FontSize',42,'TickLabelInterpreter','latex');

xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

axis equal tight
box on

if saveFigures
    save_pdf_figure(fig2,'vdp_optimal_measure.pdf');
end

fprintf('\nFinished Van der Pol moment/SOS time-average example.\n');

%% Local functions

function [objval, yval, basis] = solve_moment_sdp_vdp(d, R, opts, boundType)
% Moment relaxation for invariant probability measures.
%
% boundType = +1 maximizes Phi.
% boundType = -1 minimizes Phi.

    basis = monomial_basis_2d(d);
    N = size(basis,1);

    y = sdpvar(N,1);

    idx = make_index_map(basis);

    constraints = [];

    %% Probability normalization

    constraints = [constraints, y(get_idx(idx,0,0)) == 1];

    %% Moment matrix positivity

    rMoment = floor(d/2);
    basisMoment = monomial_basis_2d(rMoment);

    M = build_moment_matrix(y, idx, basisMoment, [0 0]);
    constraints = [constraints, M >= 0];

    %% Support on the unit disk

    rLocalizing = floor((d - 2)/2);

    if rLocalizing >= 0
        basisLocalizing = monomial_basis_2d(rLocalizing);
        Mg = build_disk_localizing_matrix(y, idx, basisLocalizing);
        constraints = [constraints, Mg >= 0];
    end

    %% Weak invariance constraints
    %
    %     int grad(p) . f dmu = 0
    %
    % for all test polynomials p of degree at most d - 2.

    degf = 3;
    testDegree = d + 1 - degf;

    if testDegree >= 1

        basisTest = monomial_basis_2d(testDegree);

        for k = 1:size(basisTest,1)

            a = basisTest(k,1);
            b = basisTest(k,2);

            if a == 0 && b == 0
                continue;
            end

            expr = 0;

            % Scaled Van der Pol dynamics:
            %
            %   X' = Y,
            %   Y' = (1 - R^2 X^2)Y - X.

            if a > 0
                expr = expr + a*moment_expr(y,idx,a-1,b+1);
            end

            if b > 0
                expr = expr ...
                    + b*moment_expr(y,idx,a,b) ...
                    - b*R^2*moment_expr(y,idx,a+2,b) ...
                    - b*moment_expr(y,idx,a+1,b-1);
            end

            constraints = [constraints, expr == 0]; %#ok<AGROW>
        end
    end

    %% Observable moment

    observableMoment = moment_expr(y,idx,2,0) ...
                     + moment_expr(y,idx,0,2);

    if boundType == +1
        objective = -observableMoment;
    elseif boundType == -1
        objective = observableMoment;
    else
        error('boundType must be +1 or -1.');
    end

    diagnostics = optimize(constraints, objective, opts);

    if diagnostics.problem ~= 0
        warning('Moment SDP diagnostic: %s', diagnostics.info);
    end

    objval = value(observableMoment);
    yval = value(y);
end

function [bound, defectNum, vars] = solve_sos_bound_vdp(d, R, opts, boundType)
% SOS auxiliary-function relaxation dual to the moment SDP.
%
% boundType = +1 computes an upper bound C satisfying
%
%     C - Phi + grad(V).f >= 0.
%
% boundType = -1 computes a lower bound C satisfying
%
%     Phi - C + grad(V).f >= 0.

    sdpvar x y C

    vars = [x; y];

    f1 = y;
    f2 = (1 - R^2*x^2)*y - x;

    Phi = x^2 + y^2;
    g = 1 - x^2 - y^2;

    degf = 3;
    Vdeg = d + 1 - degf;

    if Vdeg < 0
        error('Degree d = %d is too small for the vector field degree.', d);
    end

    [V,cV] = polynomial(vars,Vdeg);

    dVf = jacobian(V,vars)*[f1; f2];

    if boundType == +1

        defect = C - Phi + dVf;
        objective = C;

    elseif boundType == -1

        defect = Phi - C + dVf;
        objective = -C;

    else

        error('boundType must be +1 or -1.');
    end

    multiplierDegree = d - 2;

    if multiplierDegree < 0
        error('Degree d = %d is too small for the disk multiplier.', d);
    end

    [s,cs] = polynomial(vars,multiplierDegree);

    constraints = [
        sos(s), ...
        sos(defect - s*g)
        ];

    params = [cV; cs; C];

    diagnostics = solvesos(constraints, objective, opts, params);

    if diagnostics.problem ~= 0
        warning('SOS diagnostic: %s', diagnostics.info);
    end

    bound = value(C);
    defectNum = clean(replace(defect, params, value(params)), 1e-10);
end

function basis = monomial_basis_2d(d)

    basis = [];

    for total = 0:d
        for a = 0:total
            b = total - a;
            basis = [basis; a b]; %#ok<AGROW>
        end
    end
end

function idx = make_index_map(basis)

    idx = containers.Map();

    for k = 1:size(basis,1)
        key = sprintf('%d_%d', basis(k,1), basis(k,2));
        idx(key) = k;
    end
end

function k = get_idx(idx, a, b)

    key = sprintf('%d_%d', a, b);

    if ~isKey(idx,key)
        error('Moment y_(%d,%d) is outside the truncation.', a, b);
    end

    k = idx(key);
end

function expr = moment_expr(y, idx, a, b)

    expr = y(get_idx(idx,a,b));
end

function M = build_moment_matrix(y, idx, basisMoment, shift)

    n = size(basisMoment,1);
    M = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisMoment(i,1) + basisMoment(j,1) + shift(1);
            b = basisMoment(i,2) + basisMoment(j,2) + shift(2);

            M(i,j) = moment_expr(y, idx, a, b);
            M(j,i) = M(i,j);
        end
    end
end

function Mg = build_disk_localizing_matrix(y, idx, basisLocalizing)

    n = size(basisLocalizing,1);
    Mg = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisLocalizing(i,1) + basisLocalizing(j,1);
            b = basisLocalizing(i,2) + basisLocalizing(j,2);

            val = moment_expr(y,idx,a,b) ...
                - moment_expr(y,idx,a+2,b) ...
                - moment_expr(y,idx,a,b+2);

            Mg(i,j) = val;
            Mg(j,i) = val;
        end
    end
end

function [coeff, basis] = reconstruct_density_from_moments(yval, basisY, d)
% Reconstruct a polynomial density rho_d on the unit disk by matching the
% extremal moment sequence in an L^2 disk basis.

    basis = monomial_basis_2d(d);
    N = size(basis,1);

    massMatrix = zeros(N,N);
    rhs = zeros(N,1);

    idxY = make_index_map(basisY);

    for i = 1:N

        a = basis(i,1);
        b = basis(i,2);

        rhs(i) = yval(get_idx(idxY,a,b));

        for j = 1:N

            aa = basis(i,1) + basis(j,1);
            bb = basis(i,2) + basis(j,2);

            massMatrix(i,j) = disk_moment(aa,bb);
        end
    end

    coeff = massMatrix\rhs;
end

function val = disk_moment(a,b)
% Integral over the unit disk of X^a Y^b dX dY.

    if mod(a,2) == 1 || mod(b,2) == 1
        val = 0;
        return;
    end

    m = a/2;
    n = b/2;

    angular = 2*gamma(m + 0.5)*gamma(n + 0.5) ...
        / gamma(m + n + 1);

    radial = 1/(2*m + 2*n + 2);

    val = angular*radial;
end

function val = eval_poly_from_coeffs(coeff, basis, x, y)

    val = 0;

    for k = 1:length(coeff)
        val = val + coeff(k)*x^basis(k,1)*y^basis(k,2);
    end
end

function save_pdf_figure(fig,fileName)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n',fileName);
end