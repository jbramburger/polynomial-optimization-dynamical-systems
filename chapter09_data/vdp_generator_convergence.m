%% vdp_generator_convergence.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Data-Driven Auxiliary Function Methods
% Example: Pointwise convergence of learned generator approximations
%
% This script compares pointwise generator approximations for the Van der
% Pol oscillator
%
%     x' = y,
%     y' = mu(1-x^2)y - x.
%
% A trajectory on the limit cycle is sampled at several time steps tau.
% From the resulting data pairs, the script builds EDMD and gEDMD
% approximations of the infinitesimal generator and evaluates the pointwise
% error for the test observable
%
%     psi(x,y) = xy.
%
% Pressing Run should reproduce the numerical output and save:
%
%     vdp_generator_convergence_tau1.pdf
%     vdp_generator_convergence_tau2.pdf
%     vdp_generator_convergence_tau3.pdf
%     vdp_generator_convergence_tau4.pdf
%
% Requirements:
%   - MATLAB
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% User options

mu = 1;

degreePsi = 2;
degreeTheta = 4;

tauList = [0.1 0.01 0.001 0.0001];
Tdata = 600;

gridSize = 180;
xLimits = [-3 3];
yLimits = [-3 3];

observablePower = [1 1];

saveFigures = true;

%% Book plotting style

S.blue  = [0, 92, 175]/255;
S.red   = [200, 50, 50]/255;
S.black = [0, 0, 0];

%% Build dictionaries

powersPsi = monomial_powers_2d(degreePsi);
powersTheta = monomial_powers_2d(degreeTheta);

numPsi = size(powersPsi,1);
numTheta = size(powersTheta,1);

embeddingMatrix = zeros(numPsi,numTheta);

for i = 1:numPsi

    idx = find( ...
        powersTheta(:,1) == powersPsi(i,1) & ...
        powersTheta(:,2) == powersPsi(i,2), 1);

    embeddingMatrix(i,idx) = 1;
end

idxObservable = find( ...
    powersPsi(:,1) == observablePower(1) & ...
    powersPsi(:,2) == observablePower(2), 1);

if isempty(idxObservable)
    error('Requested observable is not contained in Psi dictionary.');
end

%% Get point on Van der Pol limit cycle

fprintf('\n============================================================\n');
fprintf('Van der Pol generator convergence\n');
fprintf('============================================================\n\n');

fprintf('Computing point on the Van der Pol limit cycle...\n');

x0 = [2; 0];

odeopts = odeset('RelTol',1e-11,'AbsTol',1e-13);

[~,Xtransient] = ode45(@(t,x)vdp_rhs(t,x,mu), ...
    [0 300], x0, odeopts);

xLimitCycle = Xtransient(end,:).';

fprintf('Limit-cycle point: [%.8f, %.8f]\n', ...
    xLimitCycle(1), xLimitCycle(2));

%% Evaluation grid

xv = linspace(xLimits(1),xLimits(2),gridSize);
yv = linspace(yLimits(1),yLimits(2),gridSize);

[Xg,Yg] = meshgrid(xv,yv);

Grid = [Xg(:), Yg(:)];

ThetaGrid = eval_monomials(Grid,powersTheta);

exactGenerator = exact_generator_monomial(Grid,observablePower,mu);
exactGenerator = reshape(exactGenerator,size(Xg));

%% Compute errors for all tau values

errorCells = cell(numel(tauList),1);
maxErrors = nan(size(tauList));

globalMin = inf;
globalMax = -inf;

fprintf('\nComputing learned generator errors...\n');

for q = 1:numel(tauList)

    tau = tauList(q);
    numPairs = round(Tdata/tau);

    fprintf('tau = %.4g, number of pairs = %d\n', tau, numPairs);

    [Xdata,Ydata] = generate_pairs(xLimitCycle,tau,numPairs,mu);

    [generatorEDMD,generatorGEDMD] = compute_generators( ...
        Xdata,Ydata,powersPsi,powersTheta,embeddingMatrix,tau); %#ok<ASGLU>

    approximation = generatorGEDMD(idxObservable,:)*ThetaGrid;
    approximation = reshape(approximation,size(Xg));

    Err = abs(approximation - exactGenerator);

    logErr = log10(Err + 1e-14);

    errorCells{q} = logErr;
    maxErrors(q) = max(Err(:));

    globalMin = min(globalMin,min(logErr(:)));
    globalMax = max(globalMax,max(logErr(:)));

    fprintf('  max pointwise error = %.4e\n', maxErrors(q));
end

%% Print convergence diagnostics

Results = table(tauList(:),maxErrors(:), ...
    'VariableNames',{'tau','max_pointwise_error'});

disp(Results);

%% Combined diagnostic plot

figCombined = figure;
set(figCombined,'Color','w','Units','centimeters','Position',[2 2 14 10]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for q = 1:numel(tauList)

    nexttile;

    surf(Xg,Yg,errorCells{q},'EdgeColor','none');

    view(2);
    axis equal tight;

    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');

    title(sprintf('$\\tau = %.4g$',tauList(q)), ...
        'Interpreter','latex');

    clim([globalMin globalMax]);
    colorbar;

    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
end

colormap(parula);

sgtitle( ...
    'Pointwise generator error: $\log_{10}|\widehat{\mathcal L}(xy)-\mathcal L(xy)|$', ...
    'Interpreter','latex');

if saveFigures
    exportgraphics(figCombined, ...
        'vdp_generator_convergence_combined.pdf', ...
        'ContentType','vector');
end

%% Book figures: one panel per tau value

for q = 1:numel(tauList)

    fig = figure;
    set(fig,'Color','w','Units','centimeters','Position',[2 2 11 9]);

    surf(Xg,Yg,errorCells{q},'EdgeColor','none');

    view(2);
    axis tight;
    shading interp
    box on
    grid on

    clim([globalMin globalMax]);
    colormap(parula);
    colorbar;

    set(gca,'FontSize',36,'TickLabelInterpreter','latex','LineWidth',1.0);

    xlabel('$x$','Interpreter','latex','FontSize',36);
    ylabel('$y$','Interpreter','latex','FontSize',36);

    if saveFigures
        fileName = sprintf('vdp_generator_convergence_tau%d.pdf',q);
        save_pdf_figure(fig,fileName);
    end
end

fprintf('\nFinished Van der Pol generator-convergence example.\n');

%% Local functions

function dx = vdp_rhs(~,x,mu)

    dx = [
        x(2);
        mu*(1 - x(1)^2)*x(2) - x(1)
    ];
end

function powers = monomial_powers_2d(d)

    powers = [];

    for total = 0:d
        for a = total:-1:0
            b = total - a;
            powers = [powers; a b]; %#ok<AGROW>
        end
    end
end

function Phi = eval_monomials(X,powers)

    numSnapshots = size(X,1);
    numPowers = size(powers,1);

    Phi = zeros(numPowers,numSnapshots);

    x = X(:,1).';
    y = X(:,2).';

    for j = 1:numPowers
        Phi(j,:) = x.^powers(j,1).*y.^powers(j,2);
    end
end

function Lpsi = exact_generator_monomial(X,power,mu)

    x = X(:,1).';
    y = X(:,2).';

    a = power(1);
    b = power(2);

    f1 = y;
    f2 = mu*(1 - x.^2).*y - x;

    term1 = zeros(size(x));
    term2 = zeros(size(x));

    if a > 0
        term1 = a*x.^(a-1).*y.^b.*f1;
    end

    if b > 0
        term2 = b*x.^a.*y.^(b-1).*f2;
    end

    Lpsi = term1 + term2;
end

function [Xdata,Ydata] = generate_pairs(x0,tau,numPairs,mu)

    tspan = 0:tau:numPairs*tau;

    odeopts = odeset('RelTol',1e-11,'AbsTol',1e-13);

    [~,Xsol] = ode45(@(t,x)vdp_rhs(t,x,mu),tspan,x0,odeopts);

    Xdata = Xsol(1:numPairs,:);
    Ydata = Xsol(2:numPairs+1,:);
end

function [generatorEDMD,generatorGEDMD] = compute_generators( ...
    Xdata,Ydata,powersPsi,powersTheta,embeddingMatrix,tau)

    thetaX = eval_monomials(Xdata,powersTheta);

    psiY = eval_monomials(Ydata,powersPsi);
    psiX = eval_monomials(Xdata,powersPsi);

    koopmanEDMD = psiY*pinv(thetaX);

    generatorEDMD = (koopmanEDMD - embeddingMatrix)/tau;

    finiteDifferenceData = (psiY - psiX)/tau;

    generatorGEDMD = finiteDifferenceData*pinv(thetaX);
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