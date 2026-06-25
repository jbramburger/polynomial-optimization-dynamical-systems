%% kuramoto_locked_branches.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Locked branches in finite Kuramoto rings
%
% Computes first locked-state thresholds for n = 3,...,10 and produces a
% bifurcation diagram of locked branches for one selected ring size.
%
% Output:
%   kuramoto_locked_branches_n10.pdf
%   kuramoto_locked_branches_n10_data.mat
%
% Requirements:
%   - MATLAB
%   - Optimization Toolbox

clear; clc; close all;

%% User options

nList = 3:10;
nDiag = 10;

numK = 260;
numRandomStarts = 250;

observableName = 'R2';

saveFigure = true;
saveData = true;

%% Plot style

S.blue  = [0, 92, 175]/255;
S.red   = [200, 50, 50]/255;
S.black = [0, 0, 0];
S.gray  = [0.40, 0.40, 0.40];

%% First locked-state thresholds

fprintf('\n============================================================\n');
fprintf('Kuramoto finite-ring locked branches\n');
fprintf('============================================================\n\n');

fprintf('First locked-state thresholds:\n');
fprintf(' n        K_first          stable?       max transverse eigenvalue\n');
fprintf('------------------------------------------------------------------\n');

Kfirst = nan(size(nList));
isStableFirst = false(size(nList));
maxTransEigFirst = nan(size(nList));

for idx = 1:numel(nList)

    n = nList(idx);

    omega = linspace(-1,1,n)';
    omega = omega - mean(omega);

    A = KuramotoRing.nearestNeighbourRingAdjacency(n);

    [Kc,thetaC] = KuramotoRing.firstLockedCouplingByBranchSearch(n,omega);

    if isempty(thetaC)
        fprintf('%2d        not found\n',n);
        continue
    end

    [isStable,lam] = KuramotoRing.lockedStateStability(thetaC,A,Kc);

    Kfirst(idx) = Kc;
    isStableFirst(idx) = isStable;
    maxTransEigFirst(idx) = lam(2);

    fprintf('%2d   %14.8f        %5d        %+ .6e\n', ...
        n,Kc,isStable,lam(2));
end

%% Bifurcation diagram

omega = linspace(-1,1,nDiag)';
omega = omega - mean(omega);

A = KuramotoRing.nearestNeighbourRingAdjacency(nDiag);

idxDiag = find(nList == nDiag,1);

if ~isempty(idxDiag) && ~isnan(Kfirst(idxDiag))
    KcritDiag = Kfirst(idxDiag);
else
    [KcritDiag,~] = KuramotoRing.firstLockedCouplingByBranchSearch(nDiag,omega);
end

Kmin = max(0,KcritDiag - 1.0);
Kmax = KcritDiag + 5.0;
Kgrid = linspace(Kmin,Kmax,numK);

fprintf('\nComputing branch cloud for n = %d...\n',nDiag);

branchData = KuramotoRing.computeLockedBranchCloud( ...
    nDiag,omega,A,Kgrid,numRandomStarts,1e-9,1e-5,observableName);

%% Plot

fig = figure;
set(fig,'Color','w','Units','centimeters','Position',[2 2 13 6]);

hold on
box on

stableIdx = [branchData.isStable] == true;
unstableIdx = [branchData.isStable] == false;

plot([branchData(unstableIdx).K], [branchData(unstableIdx).Y], '.', ...
    'Color',S.gray, ...
    'MarkerSize',10);

plot([branchData(stableIdx).K], [branchData(stableIdx).Y], '.', ...
    'Color',S.blue, ...
    'MarkerSize',15);

xline(KcritDiag,'--','Color',S.red,'LineWidth',3);

xlabel('$K$','Interpreter','latex','FontSize',24,'FontWeight','bold');

switch observableName
    case 'edgeCoherence'
        ylabel('$\Phi_{\rm edge}$','Interpreter','latex','FontSize',24,'FontWeight','bold');
        ylim([-0.05 1.05]);

    case 'R2'
        ylabel('order parameter','Interpreter','latex','FontSize',24);
        ylim([-0.05 1.05]);

    case 'maxEdgeDiff'
        ylabel('$\max_i |\theta_{i+1}-\theta_i|$', ...
            'Interpreter','latex','FontSize',24,'FontWeight','bold');
end

set(gca,'FontSize',24,'TickLabelInterpreter','latex','Layer','top');
grid on
xlim([min(Kgrid) max(Kgrid)]);

if saveFigure
    figName = sprintf('kuramoto_locked_branches_n%d.pdf',nDiag);
    KuramotoRing.savePDF(fig,figName);
end

if saveData
    dataName = sprintf('kuramoto_locked_branches_n%d_data.mat',nDiag);

    save(dataName, ...
        'nDiag','omega','A','Kgrid','KcritDiag','branchData', ...
        'Kfirst','nList','isStableFirst','maxTransEigFirst');

    fprintf('Saved data: %s\n',dataName);
end

fprintf('\nFinished Kuramoto locked-branch example.\n');