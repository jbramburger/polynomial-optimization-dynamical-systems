%% kuramoto_locking_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Auxiliary-function bounds for frequency locking in the Kuramoto model
%
% Computes SOS upper bounds on the locking residual for a nearest-neighbour
% Kuramoto ring with equally spaced frequencies.
%
% Output:
%   kuramoto_locking_bounds.mat
%   kuramoto_locking_bounds.pdf
%
% Requirements:
%   - MATLAB
%   - YALMIP
%   - MOSEK

clear; clc; close all;
yalmip('clear');

%% User options

n = 3;
Kgrid = 0:0.01:2;
degList = 4;

useReflectionSymmetry = false;
computeLowerBounds = false;

solverName = 'mosek';
saveFigure = true;
saveData = true;

%% Plot style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.orange = [230, 120, 20]/255;
S.purple = [170, 90, 160]/255;
S.cyan   = [0, 170, 200]/255;
S.black  = [0, 0, 0];

plotColors = {S.blue,S.red,S.green,S.purple,S.orange,S.cyan};

%% Model

omega = linspace(-1,1,n)';
omega = omega - mean(omega);

A = KuramotoRing.nearestNeighbourRingAdjacency(n);

baseCliques = KuramotoRing.buildEdgeDifferenceBaseCliques(n);
chordalCliques = KuramotoRing.chordalExtensionCliques(n,baseCliques);

fprintf('\n============================================================\n');
fprintf('Kuramoto locking-residual bounds\n');
fprintf('============================================================\n\n');

fprintf('n = %d\n',n);
fprintf('K range: %.4f to %.4f\n',min(Kgrid),max(Kgrid));
fprintf('degrees: ');
fprintf('%d ',degList);
fprintf('\n\n');

fprintf('Chordal sparse SOS cliques:\n');
for a = 1:numel(chordalCliques)
    fprintf('  clique %2d: ',a);
    fprintf('%d ',chordalCliques{a});
    fprintf('\n');
end

%% Solver settings

opts = sdpsettings('solver',solverName,'verbose',0,'sos.model',1);
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-8;
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-8;
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-8;

%% SOS sweep

Upper = nan(numel(degList),numel(Kgrid));
Lower = nan(numel(degList),numel(Kgrid));

for id = 1:numel(degList)

    d = degList(id);
    degV = max(0,d-2);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Polynomial degree d = %d, auxiliary degree <= %d\n',d,degV);
    fprintf('------------------------------------------------------------\n');

    for ik = 1:numel(Kgrid)

        K = Kgrid(ik);
        fprintf('K = %.4f ... ',K);

        %% Upper bound

        yalmip('clear');

        sv = sdpvar(n,1);
        cv = sdpvar(n,1);
        x = [sv; cv];
        h = sv.^2 + cv.^2 - 1;

        Phi = KuramotoRing.lockingResidualPolynomial(sv,cv,omega,A,K);
        f = KuramotoRing.polynomialVectorField(sv,cv,omega,A,K);

        [V,coeffV] = KuramotoRing.phaseDifferenceAuxiliary( ...
            sv,cv,degV,useReflectionSymmetry);

        C = sdpvar(1,1);
        Dupper = C - Phi - jacobian(V,x)*f;

        [Fup,coeffCertUp] = KuramotoRing.chordalTorusSOSConstraint( ...
            Dupper,h,sv,cv,d,chordalCliques);

        Fup = [Fup, C <= 1e-7]; 

        sol = solvesos(Fup,C,opts,[coeffV; coeffCertUp; C]);

        if sol.problem == 0 || sol.problem == 4
            Upper(id,ik) = value(C);
        else
            warning('Upper-bound solve failed at K=%.4f: %s',K,sol.info);
        end

        %% Optional lower bound

        if computeLowerBounds

            yalmip('clear');

            sv = sdpvar(n,1);
            cv = sdpvar(n,1);
            x = [sv; cv];
            h = sv.^2 + cv.^2 - 1;

            Phi = KuramotoRing.lockingResidualPolynomial(sv,cv,omega,A,K);
            f = KuramotoRing.polynomialVectorField(sv,cv,omega,A,K);

            [V,coeffV] = KuramotoRing.phaseDifferenceAuxiliary( ...
                sv,cv,degV,useReflectionSymmetry);

            ell = sdpvar(1,1);
            Dlower = Phi - ell + jacobian(V,x)*f;

            [Flo,coeffCertLo] = KuramotoRing.chordalTorusSOSConstraint( ...
                Dlower,h,sv,cv,d,chordalCliques);

            Flo = [Flo, ell <= 1e-7]; 

            sol = solvesos(Flo,-ell,opts,[coeffV; coeffCertLo; ell]);

            if sol.problem == 0 || sol.problem == 4
                Lower(id,ik) = value(ell);
            else
                warning('Lower-bound solve failed at K=%.4f: %s',K,sol.info);
            end
        end

        fprintf('C = %.4g\n',Upper(id,ik));
    end
end

%% Save data

if saveData
    save('kuramoto_locking_bounds.mat', ...
        'n','omega','A','Kgrid','degList','useReflectionSymmetry', ...
        'baseCliques','chordalCliques','Upper','Lower');
    fprintf('\nSaved data: kuramoto_locking_bounds.mat\n');
end

%% Plot

fig = figure;
set(fig,'Color','w','Units','centimeters','Position',[2 2 13 8]);

hold on
box on

for id = 1:numel(degList)

    col = plotColors{min(id,numel(plotColors))};

    plot(Kgrid,Upper(id,:),'-', ...
        'Color',col, ...
        'LineWidth',4);

    if computeLowerBounds
        plot(Kgrid,Lower(id,:),'--', ...
            'Color',col, ...
            'LineWidth',4);
    end
end

yline(0,'k:','LineWidth',2);

xlabel('$K$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('upper bound','Interpreter','latex','FontSize',24);

set(gca,'FontSize',24,'TickLabelInterpreter','latex','Layer','top');
grid on
xlim([min(Kgrid) max(Kgrid)]);

legendEntries = cell(1,numel(degList));
for id = 1:numel(degList)
    legendEntries{id} = sprintf('$d=%d$',degList(id));
end

legend(legendEntries, ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'Location','east');

if saveFigure
    KuramotoRing.savePDF(fig,'kuramoto_locking_bounds.pdf');
end

fprintf('\nFinished Kuramoto locking-bound example.\n');