%% shear_flow_edmd_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Data-Driven Auxiliary Function Methods
% Example: EDMD-SOS bounds for a projected shear-flow model
%
% This script computes data-driven auxiliary-function bounds for the
% nine-dimensional Moehlis shear-flow model.
%
% For each Reynolds number, the script:
%
%   1. simulates the full shear-flow model from prescribed initial data,
%   2. projects the shifted trajectories onto a low-dimensional random
%      subspace,
%   3. learns a polynomial generator approximation using EDMD, and
%   4. computes an SOS upper bound on the long-time average of ||x||^2
%      in the projected coordinates.
%
% Pressing Run computes the bounds for one auxiliary-function degree and
% saves a MAT-file of the form
%
%     shear_flow_edmd_bounds_3modes_degV5.mat
%
% Requirements:
%   - MATLAB
%   - YALMIP
%   - MOSEK
%
% Data files:
%   - ShearICs.mat
%   - random_projection_matrix.mat
%
% -------------------------------------------------------------------------

clear; clc; close all;
format long;
yalmip('clear');

%% User options

dt = 0.01;
tGrid = 0:dt:50;
ReList = 10:1:100;

numICs = 20;
numModes = 3;

degV = 5;
degF = degV + 1;

icFile = 'ShearICs.mat';
projectionFile = 'random_projection_matrix.mat';

useRowScaling = false;

rankMode = 'tol';
pinvTol = 1e-1;

if degV == 13
    pinvTol = 5e-3;
end

minRank = [];
maxRank = [];

if degV == 13
    coeffDamping = 1e-14;
else
    coeffDamping = 1e-12;
end

solver_name = 'mosek';
verbose_solver = 1;

saveFile = sprintf('shear_flow_edmd_bounds_%dmodes_degV%d.mat',numModes,degV);

%% Load data

S = load(icFile);

if isfield(S,'shearICs')
    shearICs = S.shearICs;
else
    error('Could not find variable shearICs in %s.',icFile);
end

S = load(projectionFile);

if isfield(S,'Q')
    Q = S.Q;
elseif isfield(S,'O')
    Q = S.O;
else
    error('Could not find variable Q or O in %s.',projectionFile);
end

if size(Q,2) < numModes
    error('Projection matrix has only %d columns, but numModes = %d.', ...
        size(Q,2),numModes);
end

if size(shearICs,2) < numICs
    error('shearICs has only %d initial conditions, but numICs = %d.', ...
        size(shearICs,2),numICs);
end

projectionMatrix = Q(:,1:numModes);

%% Monomial powers

generatorPowers = monomial_powers_total_degree(numModes,degF);
generatorPowers = generatorPowers(sum(generatorPowers,2) >= 1,:);

auxiliaryPowers = generatorPowers(sum(generatorPowers,2) <= degV,:);

numGeneratorMonomials = size(generatorPowers,1);
numAuxiliaryMonomials = size(auxiliaryPowers,1);

if isempty(minRank)
    minRank = numAuxiliaryMonomials;
end

if isempty(maxRank)
    maxRank = numGeneratorMonomials;
else
    maxRank = min(maxRank,numGeneratorMonomials);
end

fprintf('\n============================================================\n');
fprintf('Shear-flow EDMD-SOS bounds\n');
fprintf('============================================================\n\n');

fprintf('Number of projected modes:       %d\n',numModes);
fprintf('Auxiliary degree degV:           %d\n',degV);
fprintf('Generator dictionary degree:     %d\n',degF);
fprintf('Number of generator monomials:   %d\n',numGeneratorMonomials);
fprintf('Number of auxiliary monomials:   %d\n',numAuxiliaryMonomials);
fprintf('Rank mode:                       %s\n',rankMode);
fprintf('pinvTol:                         %.3e\n',pinvTol);
fprintf('minRank = %d, maxRank = %d\n\n',minRank,maxRank);

%% Symbolic objects

x = sdpvar(numModes,1);
generatorMonomials = monomial_vector(x,generatorPowers);

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',2);

%% Storage

upperBounds = nan(size(ReList));
diagnostics = cell(size(ReList));

retainedRank = nan(size(ReList));
conditionNumber = nan(size(ReList));
energyKept = nan(size(ReList));
tolEffective = nan(size(ReList));
coefficientNorm = nan(size(ReList));
coefficientMax = nan(size(ReList));

%% Main loop

for jj = 1:numel(ReList)

    Re = ReList(jj);

    fprintf('\n============================================================\n');
    fprintf('Re = %g  (%d of %d)\n',Re,jj,numel(ReList));
    fprintf('============================================================\n');

    %% Simulate full 9D model

    [B,L,N] = shear_flow_model(Re);
    rhs = @(~,a) B + L*a + N*reshape(a*a.',[],1);

    numSteps = numel(tGrid);
    numPairsPerIC = numSteps - 1;
    numPairs = numICs*numPairsPerIC;

    Xproj = zeros(numPairs,numModes);
    Yproj = zeros(numPairs,numModes);

    rowStart = 1;

    for ii = 1:numICs

        a0 = shearICs(:,ii);

        [~,sol] = ode45(rhs,tGrid,a0);

        % Shift laminar state to the origin.
        sol = sol - [1, zeros(1,8)];

        Xfull = sol(1:end-1,:).';
        Yfull = sol(2:end,:).';

        rows = rowStart:(rowStart+numPairsPerIC-1);

        Xproj(rows,:) = (projectionMatrix.'*Xfull).';
        Yproj(rows,:) = (projectionMatrix.'*Yfull).';

        rowStart = rowStart + numPairsPerIC;
    end

    %% Build EDMD matrices

    dictionaryX = evaluate_monomials(Xproj,generatorPowers).';
    dictionaryY = evaluate_monomials(Yproj,auxiliaryPowers).';

    %% EDMD regression

    [koopmanMatrix,regInfo] = edmd_regression( ...
        dictionaryY,dictionaryX, ...
        'useRowScaling',useRowScaling, ...
        'rankMode',rankMode, ...
        'pinvTol',pinvTol, ...
        'minRank',minRank, ...
        'maxRank',maxRank);

    retainedRank(jj) = regInfo.rank;
    conditionNumber(jj) = regInfo.cond;
    energyKept(jj) = regInfo.energyKept;
    tolEffective(jj) = regInfo.tolEffective;

    fprintf('retained rank = %d / %d\n',regInfo.rank,numGeneratorMonomials);
    fprintf('condition number = %.3e\n',regInfo.cond);
    fprintf('energy kept = %.10f\n',regInfo.energyKept);
    fprintf('effective tolerance = %.3e\n',regInfo.tolEffective);

    %% Approximate generator on auxiliary monomials

    embeddingMatrix = [eye(numAuxiliaryMonomials), ...
        zeros(numAuxiliaryMonomials,numGeneratorMonomials-numAuxiliaryMonomials)];

    generatorMatrix = (koopmanMatrix - embeddingMatrix)/dt;

    %% SOS auxiliary-function upper bound

    c = sdpvar(numAuxiliaryMonomials,1);
    C = sdpvar(1,1);

    LV = c.'*(generatorMatrix*generatorMonomials);

    defect = C - LV - dot(x,x);

    constraints = [
        sos(defect), ...
        C >= 0
        ];

    sol = solvesos( ...
        constraints, ...
        C + coeffDamping*sum(c.^2), ...
        opts, ...
        [c; C]);

    diagnostics{jj} = sol;

    if sol.problem == 0 || sol.problem == 4
        upperBounds(jj) = max(0,value(C));
    else
        upperBounds(jj) = NaN;
        warning('SOS solve failed at Re = %.3g: %s',Re,sol.info);
    end

    cval = value(c);

    coefficientNorm(jj) = norm(cval);
    coefficientMax(jj) = max(abs(cval));

    fprintf('norm(c) = %.3e\n',coefficientNorm(jj));
    fprintf('max(abs(c)) = %.3e\n',coefficientMax(jj));
    fprintf('SOS status: %s\n',sol.info);
    fprintf('Upper bound: %.15g\n',upperBounds(jj));

    save(saveFile, ...
        'ReList','upperBounds','diagnostics', ...
        'retainedRank','conditionNumber','energyKept','tolEffective', ...
        'coefficientNorm','coefficientMax','coeffDamping', ...
        'degV','degF','numModes','pinvTol','useRowScaling', ...
        'rankMode','minRank','maxRank', ...
        'icFile','projectionFile');
end

fprintf('\nSaved results to %s\n',saveFile);

%% Quick diagnostic plots

fig1 = figure;
plot(ReList,upperBounds,'.-','LineWidth',2,'MarkerSize',10);
xlabel('$Re$','Interpreter','latex','FontSize',18);
ylabel('Upper bound','Interpreter','latex','FontSize',18);
grid on

fig2 = figure;
plot(ReList,retainedRank,'.-','LineWidth',2,'MarkerSize',10);
xlabel('$Re$','Interpreter','latex','FontSize',18);
ylabel('Retained EDMD rank','Interpreter','latex','FontSize',18);
grid on

fprintf('\nFinished shear-flow EDMD-SOS computation.\n');

%% Local functions

function M = evaluate_monomials(X,powers)

    numSnapshots = size(X,1);
    numMonomials = size(powers,1);

    M = ones(numSnapshots,numMonomials);

    for j = 1:numMonomials
        for k = 1:size(X,2)
            if powers(j,k) ~= 0
                M(:,j) = M(:,j).*X(:,k).^powers(j,k);
            end
        end
    end
end

function w = monomial_vector(x,powers)

    numMonomials = size(powers,1);
    w = [];

    for j = 1:numMonomials

        term = 1;

        for k = 1:numel(x)
            if powers(j,k) ~= 0
                term = term*x(k)^powers(j,k);
            end
        end

        w = [w; term]; %#ok<AGROW>
    end
end

function powers = monomial_powers_total_degree(n,d)

    powers = [];

    for total = 0:d
        powers = [powers; fixed_degree_powers(n,total)]; %#ok<AGROW>
    end
end

function powers = fixed_degree_powers(n,d)

    if n == 1
        powers = d;
        return;
    end

    powers = [];

    for k = d:-1:0
        rest = fixed_degree_powers(n-1,d-k);
        powers = [powers; [k*ones(size(rest,1),1), rest]]; %#ok<AGROW>
    end
end

function [K,info] = edmd_regression(PsiY,PsiX,varargin)

    p = inputParser;

    addParameter(p,'useRowScaling',true);
    addParameter(p,'rankMode','tol');
    addParameter(p,'pinvTol',1e-10);
    addParameter(p,'minRank',1);
    addParameter(p,'maxRank',size(PsiX,1));

    parse(p,varargin{:});

    useRowScaling = p.Results.useRowScaling;
    rankMode = p.Results.rankMode;
    pinvTol = p.Results.pinvTol;
    minRank = p.Results.minRank;
    maxRank = p.Results.maxRank;

    if useRowScaling

        sx = sqrt(mean(PsiX.^2,2));
        sy = sqrt(mean(PsiY.^2,2));

        sx(sx == 0) = 1;
        sy(sy == 0) = 1;

        PsiXs = PsiX./sx;
        PsiYs = PsiY./sy;

    else

        sx = ones(size(PsiX,1),1);
        sy = ones(size(PsiY,1),1);

        PsiXs = PsiX;
        PsiYs = PsiY;
    end

    [U,S,V] = svd(PsiXs,'econ');
    s = diag(S);

    if isempty(s)
        error('PsiX has no singular values.');
    end

    maxAllowed = min([maxRank,numel(s),size(PsiX,1)]);

    switch lower(rankMode)

        case 'tol'
            r = sum(s > pinvTol);

        otherwise
            error('Only rankMode = ''tol'' is used in this repository script.');
    end

    r = max(r,minRank);
    r = min(r,maxAllowed);

    if r < 1
        error('No singular directions retained.');
    end

    keep = false(size(s));
    keep(1:r) = true;

    Us = U(:,keep);
    Vs = V(:,keep);
    ss = s(keep);

    Ks = PsiYs*Vs*diag(1./ss)*Us.';

    K = sy.*Ks./sx.';

    info.rank = r;
    info.singularValues = s;
    info.cond = s(1)/s(r);
    info.energyKept = sum(s(1:r).^2)/sum(s.^2);

    if r < numel(s)
        info.tolEffective = 0.5*(s(r) + s(r+1));
    else
        info.tolEffective = 0;
    end
end

function [B,L,N] = shear_flow_model(Re)

    alpha = 1/2;
    beta  = pi/2;
    gamma = 1;

    k_ab  = sqrt(alpha^2 + beta^2);
    k_ag  = sqrt(alpha^2 + gamma^2);
    k_bg  = sqrt(beta^2 + gamma^2);
    k_abg = sqrt(alpha^2 + beta^2 + gamma^2);

    lambda = [
        beta^2;
        (4/3)*beta^2 + gamma^2;
        k_bg^2;
        (1/3)*(3*alpha^2 + 4*beta^2);
        k_ab^2;
        (1/3)*(3*alpha^2 + 4*beta^2 + 3*gamma^2);
        k_abg^2;
        k_abg^2;
        9*beta^2];

    L = diag(-lambda./Re);
    B = [lambda(1)./Re; zeros(8,1)];

    N = zeros(81,9);

    T = zeros(9,9);
    T(2,3) = sqrt(3/2)*(beta*gamma)/k_bg;
    T(6,8) = -sqrt(3/2)*(beta*gamma)/k_abg;
    N(:,1) = T(:);

    T = zeros(9,9);
    T(4,6) = ((5*sqrt(2))/(3*sqrt(3)))*(gamma^2/k_ag);
    T(5,7) = -gamma^2/(sqrt(6)*k_ag);
    T(5,8) = -(alpha*beta*gamma)/(sqrt(6)*k_ag*k_abg);
    T(1,3) = -sqrt(3/2)*beta*gamma/k_bg;
    T(3,9) = -sqrt(3/2)*beta*gamma/k_bg;
    N(:,2) = T(:);

    T = zeros(9,9);
    T(4,7) = (2/sqrt(6))*(alpha*beta*gamma)/(k_ag*k_bg);
    T(5,6) = (2/sqrt(6))*(alpha*beta*gamma)/(k_ag*k_bg);
    T(4,8) = (beta^2*(3*alpha^2 + gamma^2) ...
        - 3*gamma^2*(alpha^2+gamma^2))/(sqrt(6)*k_ag*k_bg*k_abg);
    N(:,3) = T(:);

    T = zeros(9,9);
    T(1,5) = -alpha/sqrt(6);
    T(5,9) = -alpha/sqrt(6);
    T(2,6) = -(10/(3*sqrt(6)))*(alpha^2/k_ag);
    T(3,7) = -sqrt(3/2)*(alpha*beta*gamma)/(k_ag*k_bg);
    T(3,8) = -sqrt(3/2)*(alpha^2*beta^2)/(k_ag*k_bg*k_abg);
    N(:,4) = T(:);

    T = zeros(9,9);
    T(1,4) = alpha/sqrt(6);
    T(4,9) = alpha/sqrt(6);
    T(3,6) = (2/sqrt(6))*(alpha*beta*gamma)/(k_ag*k_bg);
    T(2,7) = alpha^2/(sqrt(6)*k_ag);
    T(2,8) = -(alpha*beta*gamma)/(sqrt(6)*k_ag*k_abg);
    N(:,5) = T(:);

    T = zeros(9,9);
    T(1,7) = alpha/sqrt(6);
    T(1,8) = sqrt(3/2)*(beta*gamma)/k_abg;
    T(2,4) = (10/(3*sqrt(6)))*((alpha^2 - gamma^2)/k_ag);
    T(3,5) = -2*sqrt(2/3)*(alpha*beta*gamma)/(k_ag*k_bg);
    T(7,9) = alpha/sqrt(6);
    T(8,9) = sqrt(3/2)*(beta*gamma)/k_abg;
    N(:,6) = T(:);

    T = zeros(9,9);
    T(1,6) = -alpha/sqrt(6);
    T(6,9) = -alpha/sqrt(6);
    T(2,5) = (1/sqrt(6))*((gamma^2 - alpha^2)/k_ag);
    T(3,4) = (1/sqrt(6))*(alpha*beta*gamma)/(k_ag*k_bg);
    N(:,7) = T(:);

    T = zeros(9,9);
    T(2,5) = (2/sqrt(6))*(alpha*beta*gamma)/(k_ag*k_abg);
    T(3,4) = gamma^2*(3*alpha^2 - beta^2 + 3*gamma^2) ...
        /(sqrt(6)*k_ag*k_bg*k_abg);
    N(:,8) = T(:);

    T = zeros(9,9);
    T(2,3) = sqrt(3/2)*beta*gamma/k_bg;
    T(6,8) = -sqrt(3/2)*beta*gamma/k_abg;
    N(:,9) = T(:);

    N = N.';
end