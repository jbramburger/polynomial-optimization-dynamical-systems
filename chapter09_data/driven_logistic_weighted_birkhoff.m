%% driven_logistic_weighted_birkhoff.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Data-Driven Auxiliary Function Methods
% Example: Weighted Birkhoff averages for a driven logistic map
%
% This script compares ordinary and weighted Birkhoff averages for the
% quasiperiodically driven logistic map
%
%     theta_{n+1} = theta_n + omega mod 1,
%     x_{n+1}     = 3.5(1 + eps cos(2*pi*theta_n)) x_n(1-x_n).
%
% The examples compare three regimes:
%
%   eps = 0.00   periodic response,
%   eps = 0.01   quasiperiodic response,
%   eps = 0.10   chaotic response.
%
% Pressing Run should reproduce the numerical output and save:
%
%     logistic_periodic.pdf
%     logistic_quasiperiodic.pdf
%     logistic_chaotic.pdf
%
% Requirements:
%   - MATLAB
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% User options

omega = (sqrt(5) - 1)/2;

epsVals = [0.00, 0.01, 0.10];
epsTags = {'periodic','quasiperiodic','chaotic'};

Nmax = 1e6;
Nref = 1e7;
Ndiscard = 100;

Nlist = unique(round(logspace(2,log10(Nmax),1000)));

x0 = 0.25;
theta0 = 0;

saveFigures = true;

%% Book plotting style

S.blue  = [0, 92, 175]/255;
S.red   = [200, 50, 50]/255;
S.black = [0, 0, 0];

%% Run experiments

results = struct();

fprintf('\n============================================================\n');
fprintf('Weighted Birkhoff averages for the driven logistic map\n');
fprintf('============================================================\n\n');

for r = 1:numel(epsVals)

    epsVal = epsVals(r);

    fprintf('Running eps = %.3f (%s)...\n',epsVal,epsTags{r});

    [x,theta] = driven_logistic(epsVal,omega,x0,theta0,Nref+Ndiscard);

    x = x(Ndiscard+1:end);
    theta = theta(Ndiscard+1:end); %#ok<NASGU>

    observable = x;

    if epsVal <= 0.02
        wref = endpoint_flat_weights(Nref);
        Aref = sum(wref.*observable(:));
    else
        Aref = mean(observable);
    end

    obsShort = observable(1:Nmax);

    errUnweighted = zeros(size(Nlist));
    errWeighted   = zeros(size(Nlist));

    for q = 1:numel(Nlist)

        N = Nlist(q);
        hN = obsShort(1:N);

        Aunweighted = mean(hN);

        w = endpoint_flat_weights(N);
        Aweighted = sum(w.*hN(:));

        errUnweighted(q) = abs(Aunweighted - Aref);
        errWeighted(q)   = abs(Aweighted   - Aref);
    end

    guide = errUnweighted(1)*(Nlist(1)./Nlist);

    results(r).eps = epsVal;
    results(r).tag = epsTags{r};
    results(r).Nlist = Nlist;
    results(r).errUnweighted = errUnweighted;
    results(r).errWeighted = errWeighted;
    results(r).guide = guide;
    results(r).Aref = Aref;
    results(r).xplot = x(1:5000);

    fprintf('  reference average = %.16e\n',Aref);
end

%% Plot convergence figures

plot_convergence(results(1),S, ...
    [1e-16 1e-2], ...
    'logistic_periodic.pdf', ...
    saveFigures);

plot_convergence(results(2),S, ...
    [1e-16 1e-2], ...
    'logistic_quasiperiodic.pdf', ...
    saveFigures);

plot_convergence(results(3),S, ...
    [1e-7 3e-2], ...
    'logistic_chaotic.pdf', ...
    saveFigures);

fprintf('\nFinished driven logistic weighted-Birkhoff example.\n');

%% Local functions

function [x,theta] = driven_logistic(epsVal,omega,x0,theta0,N)

    x = zeros(N,1);
    theta = zeros(N,1);

    x(1) = x0;
    theta(1) = mod(theta0,1);

    for n = 1:N-1

        a = 3.5*(1 + epsVal*cos(2*pi*theta(n)));

        x(n+1) = a*x(n)*(1-x(n));
        theta(n+1) = mod(theta(n) + omega,1);

        if x(n+1) < 0
            x(n+1) = 0;
        elseif x(n+1) > 1
            x(n+1) = 1;
        end
    end
end

function w = endpoint_flat_weights(N)

    k = (1:N)';
    t = k/(N+1);

    w = exp(-1./(t.*(1-t)));
    w = w/sum(w);
end

function plot_convergence(result,S,ylims,fileName,saveFigure)

    fig = figure;
    set(fig,'Color','w','Units','centimeters','Position',[2 2 11 8]);

    errFloor = 1e-16;

    errUnweightedPlot = max(result.errUnweighted, errFloor);
    errWeightedPlot   = max(result.errWeighted,   errFloor);
    guidePlot         = max(result.guide,         errFloor);

    hold on
    box on
    grid on

    loglog(result.Nlist,errUnweightedPlot, ...
        'Color',S.blue, ...
        'LineWidth',3);

    loglog(result.Nlist,errWeightedPlot, ...
        'Color',S.red, ...
        'LineWidth',3);

    loglog(result.Nlist,guidePlot, ...
        'k--', ...
        'LineWidth',2);

    set(gca,'XScale','log','YScale','log');

    xlabel('$N$','Interpreter','latex','FontSize',24);
    ylabel('absolute error','Interpreter','latex','FontSize',24);

    ylim(ylims);

    set(gca,'FontSize',18,'TickLabelInterpreter','latex');

    legend({'unweighted','weighted','$1/N$'}, ...
        'Interpreter','latex', ...
        'Location','southwest');

    if saveFigure
        save_pdf_figure(fig,fileName);
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