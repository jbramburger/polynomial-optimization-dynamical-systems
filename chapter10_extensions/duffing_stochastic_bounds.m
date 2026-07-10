%% duffing_stochastic_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Extensions and Frontiers
% Example: Moment-SDP and SOS bounds for a stochastic Duffing oscillator
%
% This script computes moment-SDP and SOS auxiliary-function bounds on the
% stationary expectation
%
%     E[x^2 + y^2]
%
% for the stochastic Duffing/Langevin oscillator
%
%     dx = y dt,
%     dy = (x - x^3 - gamma*y) dt + sqrt(2*gamma*epsilon) dW.
%
% The computation is performed in scaled variables
%
%     x = R X,     y = R Y,
%
% on the unit disk X^2 + Y^2 <= 1. The scaled generator is
%
%     L V =
%         Y V_X
%       + (X - R^2 X^3 - gamma Y) V_Y
%       + (gamma epsilon/R^2) V_YY.
%
% Pressing Run should reproduce the numerical output and save:
%
%     duffing_stochastic_bounds.csv
%     duffing_time_average_defect.pdf
%     duffing_optimal_measure.pdf
%     duffing_exact_invariant_density.pdf
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

gammaD = 1.0;
epsilonD = 0.3;

R = 4.0;

degrees = 2:2:16;
plotDegree = 16;

gridSize = 512;

solverName = 'mosek';
verboseSolver = 0;

saveFigures = true;
saveTable = true;

%% Book plotting style

S.blue  = [0, 92, 175]/255;
S.red   = [200, 50, 50]/255;
S.black = [0, 0, 0];

%% Solver options

optsMoment = sdpsettings('solver',solverName,'verbose',verboseSolver);
optsSOS    = sdpsettings('solver',solverName,'verbose',verboseSolver);

%% Storage

numDegrees = numel(degrees);

momentLower = nan(numDegrees,1);
momentUpper = nan(numDegrees,1);
sosLower    = nan(numDegrees,1);
sosUpper    = nan(numDegrees,1);

upperMomentStore = cell(numDegrees,1);
momentBasisStore = cell(numDegrees,1);

upperDefectStore = cell(numDegrees,1);
defectVarsStore  = cell(numDegrees,1);

%% Reference values

fprintf('\n============================================================\n');
fprintf('Stochastic Duffing oscillator: moment and SOS bounds\n');
fprintf('============================================================\n\n');

fprintf('gamma      = %.6g\n', gammaD);
fprintf('epsilon    = %.6g\n', epsilonD);
fprintf('R          = %.6g\n', R);
fprintf('degrees    = ');
fprintf('%d ',degrees);
fprintf('\n\n');

exactMean = exact_stationary_mean(gammaD,epsilonD);
simulationMean = simulate_stationary_mean(gammaD,epsilonD);

fprintf('Exact full-space stationary E[x^2+y^2] = %.10f\n', exactMean);
fprintf('Euler--Maruyama estimate               = %.10f\n\n', simulationMean);

%% Compute moment and SOS hierarchies

for kk = 1:numDegrees

    d = degrees(kk);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Relaxation degree d = %d\n', d);
    fprintf('------------------------------------------------------------\n');

    try
        [momentUpper(kk), yUpper, basisY] = moment_bound( ...
            d,R,gammaD,epsilonD,optsMoment,+1);

        [momentLower(kk), ~, ~] = moment_bound( ...
            d,R,gammaD,epsilonD,optsMoment,-1);

        upperMomentStore{kk} = yUpper;
        momentBasisStore{kk} = basisY;

        fprintf('Moment lower/upper, scaled: %.10f   %.10f\n', ...
            momentLower(kk),momentUpper(kk));

    catch ME
        warning('Moment SDP failed at degree %d:\n%s',d,ME.message);
    end

    try
        [sosUpper(kk), defectUpper, defectVars] = sos_bound( ...
            d,R,gammaD,epsilonD,optsSOS,+1);

        [sosLower(kk), ~, ~] = sos_bound( ...
            d,R,gammaD,epsilonD,optsSOS,-1);

        upperDefectStore{kk} = defectUpper;
        defectVarsStore{kk} = defectVars;

        fprintf('SOS lower/upper, scaled:    %.10f   %.10f\n', ...
            sosLower(kk),sosUpper(kk));

    catch ME
        warning('SOS problem failed at degree %d:\n%s',d,ME.message);
    end
end

%% Summary table in original coordinates

scaleObservable = R^2;

Results = table( ...
    degrees(:), ...
    scaleObservable*momentLower, ...
    scaleObservable*momentUpper, ...
    scaleObservable*sosLower, ...
    scaleObservable*sosUpper, ...
    'VariableNames', ...
    {'Degree','MomentLower','MomentUpper','SOSLower','SOSUpper'});

fprintf('\n============================================================\n');
fprintf('Summary in original coordinates\n');
fprintf('============================================================\n\n');

fprintf('Exact full-space stationary E[x^2+y^2] = %.10f\n', exactMean);
fprintf('Euler--Maruyama estimate               = %.10f\n\n', simulationMean);

disp(Results);

if saveTable
    writetable(Results,'duffing_stochastic_bounds.csv');
    fprintf('Saved table: duffing_stochastic_bounds.csv\n');
end

%% Reconstruct extremal invariant density and SOS defect

plotIndex = find(degrees == plotDegree,1);

if isempty(plotIndex)
    error('Requested plotDegree = %d is not included in degrees.',plotDegree);
end

if isempty(upperMomentStore{plotIndex}) || isempty(upperDefectStore{plotIndex})
    error('Degree %d did not solve successfully, so plots cannot be generated.',plotDegree);
end

fprintf('\nReconstructing density and defect at degree d = %d...\n',plotDegree);

yExtremal = upperMomentStore{plotIndex};
basisY = momentBasisStore{plotIndex};

[densityCoeff,densityBasis] = reconstruct_density( ...
    yExtremal,basisY,plotDegree);

defectPolynomial = upperDefectStore{plotIndex};
defectVars = defectVarsStore{plotIndex};

%% Evaluate fields on grid

xx = linspace(-1,1,gridSize);
yy = linspace(-1,1,gridSize);

[Xg,Yg] = meshgrid(xx,yy);
mask = Xg.^2 + Yg.^2 <= 1;

rhoMomentGrid = nan(size(Xg));
rhoExactGrid  = nan(size(Xg));
defectGrid    = nan(size(Xg));

for q = 1:numel(Xg)

    if mask(q)

        rhoMomentGrid(q) = evaluate_polynomial( ...
            densityCoeff,densityBasis,Xg(q),Yg(q));

        defectGrid(q) = value(replace( ...
            defectPolynomial,defectVars,[Xg(q);Yg(q)]));
    end
end

rhoTmp = invariant_density_duffing(R*Xg,R*Yg,epsilonD);
rhoExactGrid(mask) = rhoTmp(mask);

massTmp = trapz(yy,trapz(xx,nan_to_zero(rhoExactGrid),2));
rhoExactGrid = rhoExactGrid/massTmp;

defectGrid(defectGrid < 1e-14) = 1e-14;

%% Figure 1: SOS defect

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 11 9]);

contourf(R*Xg,R*Yg,log(defectGrid),35,'LineColor','none');
colorbar;
colormap(parula);

set(gca,'FontSize',42,'TickLabelInterpreter','latex');

xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

axis equal tight
box on

if saveFigures
    export_pdf(fig1,'duffing_time_average_defect.pdf');
end

%% Figure 2: reconstructed moment density

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 11 9]);

contourf(R*Xg,R*Yg,rhoMomentGrid,35,'LineColor','none');
colorbar;
colormap(parula);

set(gca,'FontSize',42,'TickLabelInterpreter','latex');

xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

axis equal tight
box on

if saveFigures
    export_pdf(fig2,'duffing_optimal_measure.pdf');
end

%% Figure 3: exact invariant density

fig3 = figure;
set(fig3,'Color','w','Units','centimeters','Position',[2 2 11 9]);

contourf(R*Xg,R*Yg,rhoExactGrid,35,'LineColor','none');
colorbar;
colormap(parula);

set(gca,'FontSize',42,'TickLabelInterpreter','latex');

xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

axis equal tight
box on

if saveFigures
    export_pdf(fig3,'duffing_exact_invariant_density.pdf');
end

fprintf('\nFinished stochastic Duffing moment/SOS example.\n');

%% Local functions

function [objectiveValue,yValue,basis] = moment_bound( ...
    d,R,gammaD,epsilonD,opts,boundType)
% Moment relaxation for stationary probability measures.
%
% boundType = +1 maximizes Phi.
% boundType = -1 minimizes Phi.

    basis = monomial_basis_2d(d);
    numMoments = size(basis,1);

    y = sdpvar(numMoments,1);

    idx = make_index_map(basis);

    constraints = [];

    %% Probability normalization

    constraints = [constraints, y(get_index(idx,0,0)) == 1];

    %% Moment matrix positivity

    momentDegree = floor(d/2);
    basisMoment = monomial_basis_2d(momentDegree);

    M = build_moment_matrix(y,idx,basisMoment,[0 0]);
    constraints = [constraints, M >= 0];

    %% Unit-disk support

    localizingDegree = floor((d-2)/2);

    if localizingDegree >= 0
        basisLocalizing = monomial_basis_2d(localizingDegree);
        Mg = build_ball_localizing_matrix(y,idx,basisLocalizing);
        constraints = [constraints, Mg >= 0];
    end

    %% Stationarity constraints
    %
    %     int L p dmu = 0
    %
    % for all polynomial test functions p of degree at most d - 2.

    testDegree = d - 2;

    if testDegree >= 1

        basisTest = monomial_basis_2d(testDegree);

        for k = 1:size(basisTest,1)

            a = basisTest(k,1);
            b = basisTest(k,2);

            if a == 0 && b == 0
                continue;
            end

            expr = 0;

            % p_X * Y.
            if a > 0
                expr = expr + a*moment_expr(y,idx,a-1,b+1);
            end

            % p_Y * (X - R^2 X^3 - gamma Y).
            if b > 0
                expr = expr ...
                    + b*moment_expr(y,idx,a+1,b-1) ...
                    - b*R^2*moment_expr(y,idx,a+3,b-1) ...
                    - b*gammaD*moment_expr(y,idx,a,b);
            end

            % Diffusion term: (gamma epsilon/R^2) p_YY.
            if b >= 2
                expr = expr ...
                    + (gammaD*epsilonD/R^2)*b*(b-1) ...
                    * moment_expr(y,idx,a,b-2);
            end

            constraints = [constraints, expr == 0]; 
        end
    end

    %% Observable Phi_scaled = X^2 + Y^2

    observableMoment = moment_expr(y,idx,2,0) ...
                     + moment_expr(y,idx,0,2);

    if boundType == +1
        objective = -observableMoment;
    elseif boundType == -1
        objective = observableMoment;
    else
        error('boundType must be +1 or -1.');
    end

    diagnostics = optimize(constraints,objective,opts);

    if diagnostics.problem ~= 0
        warning('Moment SDP diagnostic: %s',diagnostics.info);
    end

    objectiveValue = value(observableMoment);
    yValue = value(y);
end

function [bound,defectNum,vars] = sos_bound( ...
    d,R,gammaD,epsilonD,opts,boundType)
% SOS auxiliary-function relaxation.
%
% boundType = +1 computes an upper bound C satisfying
%
%     C - Phi + L V >= 0.
%
% boundType = -1 computes a lower bound C satisfying
%
%     Phi - C + L V >= 0.

    X = sdpvar(1,1);
    Y = sdpvar(1,1);

    vars = [X; Y];

    Phi = X^2 + Y^2;
    g = 1 - X^2 - Y^2;

    Vdegree = d - 2;

    if Vdegree < 0
        error('Degree d = %d is too small.',d);
    end

    [V,cV] = polynomial(vars,Vdegree);

    C = sdpvar(1,1);

    LV = duffing_generator(V,vars,R,gammaD,epsilonD);

    if boundType == +1

        defect = C - Phi + LV;
        objective = C;

    elseif boundType == -1

        defect = Phi - C + LV;
        objective = -C;

    else

        error('boundType must be +1 or -1.');
    end

    multiplierDegree = d - 2;

    if multiplierDegree < 0
        error('Degree d = %d is too small for support multiplier.',d);
    end

    [s1,cs1] = polynomial(vars,multiplierDegree);

    constraints = [
        sos(s1), ...
        sos(defect - s1*g)
        ];

    params = [cV; cs1; C];

    diagnostics = solvesos(constraints,objective,opts,params);

    if diagnostics.problem ~= 0
        warning('SOS diagnostic: %s',diagnostics.info);
    end

    bound = value(C);
    defectNum = clean(replace(defect,params,value(params)),1e-10);
end

function LV = duffing_generator(V,vars,R,gammaD,epsilonD)

    X = vars(1);
    Y = vars(2);

    f1 = Y;
    f2 = X - R^2*X^3 - gammaD*Y;

    gradV = jacobian(V,vars);
    VYY = jacobian(jacobian(V,Y),Y);

    LV = gradV*[f1; f2] + (gammaD*epsilonD/R^2)*VYY;
end

function meanPhi = exact_stationary_mean(gammaD,epsilonD) 
% Exact full-space stationary density:
%
%     rho(x,y) proportional to exp(-(y^2/2 + U(x))/epsilon),
%
% where U(x) = x^4/4 - x^2/2. Since y is Gaussian,
%
%     E[y^2] = epsilon.

    U = @(x) 0.25*x.^4 - 0.5*x.^2;

    Zx = integral(@(x) exp(-U(x)/epsilonD), -Inf, Inf, ...
        'RelTol',1e-10, ...
        'AbsTol',1e-12);

    Ex2 = integral(@(x) x.^2.*exp(-U(x)/epsilonD), -Inf, Inf, ...
        'RelTol',1e-10, ...
        'AbsTol',1e-12)/Zx;

    Ey2 = epsilonD;

    meanPhi = Ex2 + Ey2;
end

function rho = invariant_density_duffing(X,Y,epsilonD)

    U = 0.25*X.^4 - 0.5*X.^2;
    H = 0.5*Y.^2 + U;

    rho = exp(-H/epsilonD);
end

function meanPhi = simulate_stationary_mean(gammaD,epsilonD)
% Long Euler--Maruyama estimate. This is only a sanity check and is not
% used in the optimization.

    rng(1);

    dt = 1e-3;
    Tburn = 200;
    Tavg = 1000;

    numBurn = round(Tburn/dt);
    numAvg = round(Tavg/dt);

    x = 1.0;
    y = 0.0;

    sigmaNoise = sqrt(2*gammaD*epsilonD);

    for k = 1:numBurn 

        xOld = x;
        yOld = y;

        x = xOld + dt*yOld;

        y = yOld ...
            + dt*(xOld - xOld^3 - gammaD*yOld) ...
            + sigmaNoise*sqrt(dt)*randn;
    end

    acc = 0;

    for k = 1:numAvg 

        xOld = x;
        yOld = y;

        x = xOld + dt*yOld;

        y = yOld ...
            + dt*(xOld - xOld^3 - gammaD*yOld) ...
            + sigmaNoise*sqrt(dt)*randn;

        acc = acc + (x^2 + y^2);
    end

    meanPhi = acc/numAvg;
end

function basis = monomial_basis_2d(d)

    basis = [];

    for total = 0:d
        for a = 0:total
            b = total - a;
            basis = [basis; a b];
        end
    end
end

function idx = make_index_map(basis)

    idx = containers.Map();

    for k = 1:size(basis,1)
        key = sprintf('%d_%d',basis(k,1),basis(k,2));
        idx(key) = k;
    end
end

function k = get_index(idx,a,b)

    key = sprintf('%d_%d',a,b);

    if ~isKey(idx,key)
        error('Moment y_(%d,%d) is outside the truncation.',a,b);
    end

    k = idx(key);
end

function expr = moment_expr(y,idx,a,b)

    expr = y(get_index(idx,a,b));
end

function M = build_moment_matrix(y,idx,basisMoment,shift)

    n = size(basisMoment,1);

    M = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisMoment(i,1) + basisMoment(j,1) + shift(1);
            b = basisMoment(i,2) + basisMoment(j,2) + shift(2);

            M(i,j) = moment_expr(y,idx,a,b);
            M(j,i) = M(i,j);
        end
    end
end

function Mg = build_ball_localizing_matrix(y,idx,basisLocalizing)

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

function [coeff,basis] = reconstruct_density(yValue,basisY,d)
% Construct rho_d(X,Y) = coeff' z_d(X,Y) such that
%
%     int_disk z_d rho_d dX dY
%
% matches the moment sequence y.

    basis = monomial_basis_2d(d);
    numBasis = size(basis,1);

    massMatrix = zeros(numBasis,numBasis);
    rhs = zeros(numBasis,1);

    idxY = make_index_map(basisY);

    for i = 1:numBasis

        a = basis(i,1);
        b = basis(i,2);

        rhs(i) = yValue(get_index(idxY,a,b));

        for j = 1:numBasis

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

function val = evaluate_polynomial(coeff,basis,x,y)

    val = 0;

    for k = 1:length(coeff)
        val = val + coeff(k)*x^basis(k,1)*y^basis(k,2);
    end
end

function A = nan_to_zero(A)

    A(isnan(A)) = 0;
end

function export_pdf(fig,fileName)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n',fileName);
end