%% ks_background.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Background-method bounds for the Kuramoto--Sivashinsky equation
%
% This script compares background-method bounds and polynomial
% auxiliary function SOS bounds for a Galerkin truncation of the
% Kuramoto--Sivashinsky equation on the odd subspace.
%
% Pressing Run should reproduce the numerical output and save:
%
%     ks_bounds_sweep.mat
%     ks_background_L4.pdf
%     ks_background_bounds.pdf
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

N     = 16;
Mbg   = 16;
Lvals = 0:0.01:4;
Lplot = 4.0;

runBackground = true;
runDeg2       = false;
runDeg3       = false;
runDeg4       = true;

computeSteadyStates = false;

saveFigures = true;
saveData    = true;

outfile = 'ks_bounds_sweep.mat';

solver_name = 'mosek';
verbose_solver = 0;

ops = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'sos.model',1);

%% Book plotting style

S.blue  = [0, 92, 175]/255;
S.red   = [200, 50, 50]/255;
S.gray  = [0.40, 0.40, 0.40];

%% Storage

nL = numel(Lvals);

C_bg = nan(nL,1);
C_d2 = nan(nL,1);
C_d3 = nan(nL,1);
C_d4 = nan(nL,1);

a_bg_store   = nan(Mbg,nL);
eta_bg_store = nan(nL,1);

steady_E = cell(nL,1);
steady_B = cell(nL,1);

fprintf('\n============================================================\n');
fprintf('Kuramoto--Sivashinsky background-method bounds\n');
fprintf('============================================================\n\n');

fprintf('Galerkin modes N = %d\n', N);
fprintf('Background modes M = %d\n', Mbg);
fprintf('L range: %.2f to %.2f\n\n', min(Lvals), max(Lvals));

%% Sweep over L

for iL = 1:nL

    L = Lvals(iL);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('L = %.2f (%d/%d)\n', L, iL, nL);
    fprintf('------------------------------------------------------------\n');

    if L == 0
        C_bg(iL) = 0;
        C_d2(iL) = 0;
        C_d3(iL) = 0;
        C_d4(iL) = 0;
        steady_E{iL} = 0;
        steady_B{iL} = zeros(1,N);
        continue;
    end

    if runBackground
        try
            [C_bg(iL), a_bg_store(:,iL), eta_bg_store(iL)] = ...
                solve_background_KS(N,Mbg,L,ops);

            fprintf('  background method: %.8g\n', C_bg(iL));
        catch ME
            warning('Background solve failed at L=%.3f: %s', L, ME.message);
        end
    end

    if runDeg2
        try
            C_d2(iL) = solve_aux_sos_KS(N,L,2,ops);
            fprintf('  degree 2 SOS:     %.8g\n', C_d2(iL));
        catch ME
            warning('Degree 2 SOS failed at L=%.3f: %s', L, ME.message);
        end
    end

    if runDeg3
        try
            C_d3(iL) = solve_aux_sos_KS(N,L,3,ops);
            fprintf('  degree 3 SOS:     %.8g\n', C_d3(iL));
        catch ME
            warning('Degree 3 SOS failed at L=%.3f: %s', L, ME.message);
        end
    end

    if runDeg4
        try
            C_d4(iL) = solve_aux_sos_KS(N,L,4,ops);
            fprintf('  degree 4 SOS:     %.8g\n', C_d4(iL));
        catch ME
            warning('Degree 4 SOS failed at L=%.3f: %s', L, ME.message);
        end
    end

    if computeSteadyStates
        try
            [Bss, Ess] = find_steady_states_KS(N,L);
            steady_B{iL} = Bss;
            steady_E{iL} = Ess;
        catch ME
            warning('Steady-state search failed at L=%.3f: %s', L, ME.message);
        end
    end

    if saveData
        save(outfile, ...
            'Lvals','N','Mbg', ...
            'C_bg','C_d2','C_d3','C_d4', ...
            'a_bg_store','eta_bg_store', ...
            'steady_E','steady_B');
    end
end

if saveData
    fprintf('\nSaved data: %s\n', outfile);
end

%% Plot optimized background field

[~,idxPlot] = min(abs(Lvals-Lplot));
L0 = Lvals(idxPlot);

xgrid = linspace(0,2*pi*L0,1000);
tau = zeros(size(xgrid));

if all(isnan(a_bg_store(:,idxPlot)))
    warning('No background coefficients available at L = %.2f. Skipping tau plot.', L0);
else
    for j = 1:Mbg
        tau = tau + a_bg_store(j,idxPlot)*sin(j*xgrid/L0);
    end

    fig1 = figure;
    set(fig1,'Color','w','Units','centimeters','Position',[2 2 13 8]);

    plot(xgrid,tau,'Color',S.blue,'LineWidth',7);

    xlabel('$x$','Interpreter','latex','FontSize',24);
    ylabel('$\tau(x)$','Interpreter','latex','FontSize',24);

    set(gca,'FontSize',24);
    xlim([0 2*pi*L0])
    grid on
    box on

    if saveFigures
        save_pdf_figure(fig1,'ks_background_L4.pdf');
    end
end

%% Plot energy bounds and steady states

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 13 8]);

hold on
box on

if computeSteadyStates
    for iL = 1:nL
        Ess = steady_E{iL};
        if ~isempty(Ess)
            plot(Lvals(iL)*ones(size(Ess)),Ess,'.', ...
                'Color',S.gray, ...
                'MarkerSize',20);
        end
    end
end

if any(~isnan(C_d2))
    plot(Lvals,C_d2,'Color',[0.65 0.65 0.65],'LineWidth',3);
end

if any(~isnan(C_d3))
    plot(Lvals,C_d3,'Color',[0.35 0.35 0.35],'LineWidth',3);
end

if any(~isnan(C_d4))
    plot(Lvals,C_d4,'Color',S.red,'LineWidth',7);
end

if any(~isnan(C_bg))
    plot(Lvals,C_bg,'Color',S.blue,'LineWidth',7);
end

xlabel('$L$','Interpreter','latex','FontSize',24);
ylabel('$\overline{\langle u^2\rangle}$','Interpreter','latex','FontSize',24);

set(gca,'FontSize',24);
grid on

if saveFigures
    save_pdf_figure(fig2,'ks_background_bounds.pdf');
end

fprintf('\nFinished KS background-method example.\n');

%% Local functions

function [Cval, tauval, etaval, rhoval] = solve_background_KS(N,M,L,ops)

    b   = sdpvar(N,1);
    rho = sdpvar(M,1);
    eta = sdpvar(1);
    C   = sdpvar(1);

    nq = 4*(2*N + M) + 64;
    th = (0:nq-1)'*(2*pi/nq);

    S   = zeros(nq,N);
    Cx  = zeros(nq,N);
    Sxx = zeros(nq,N);

    for k = 1:N
        S(:,k)   = sin(k*th);
        Cx(:,k)  = (k/L)*cos(k*th);
        Sxx(:,k) = -(k/L)^2*sin(k*th);
    end

    Rx  = zeros(nq,M);
    Rxx = zeros(nq,M);

    for j = 1:M
        Rx(:,j)  = (j/L)*cos(j*th);
        Rxx(:,j) = -(j/L)^2*sin(j*th);
    end

    u   = S*b;
    ux  = Cx*b;
    uxx = Sxx*b;

    rho_x  = Rx*rho;
    rho_xx = Rxx*rho;

    avg = @(z) sum(z)/nq;

    Phi = avg(u.^2);

    Vdot = eta*(avg(ux.^2) - avg(uxx.^2)) ...
           - 0.5*avg(rho_x.*u.^2) ...
           - avg(rho_x.*ux) ...
           + avg(rho_xx.*uxx);

    defect = clean(C - Phi - Vdot,1e-10);

    b0 = zeros(N,1);
    const = replace(defect,b,b0);
    ell   = replace(jacobian(defect,b).',b,b0);
    Q     = 0.5*hessian(defect,b);

    G = clean([const, 0.5*ell.'; 0.5*ell, Q],1e-10);

    constraints = [
        G >= 0, ...
        C >= 0, ...
        eta >= 0
        ];

    sol = optimize(constraints,C,ops);

    if sol.problem ~= 0
        warning('Background LMI status: %s', sol.info);
    end

    Cval   = value(C);
    etaval = value(eta);
    rhoval = value(rho);

    if etaval > 1e-8
        tauval = rhoval/etaval;
    else
        tauval = zeros(M,1);
    end
end

function Cval = solve_aux_sos_KS(N,L,degV,ops)

    x = sdpvar(N,1);
    C = sdpvar(1);

    f = ks_galerkin_sdp(x,L);

    Phi = 0.5*(x.'*x);

    exps = invariant_exponents_half_shift(N,degV);
    exps = exps(sum(exps,2) > 0,:);

    mon = monomials_from_exponents(x,exps);
    cV = sdpvar(length(mon),1);

    V = cV.'*mon;

    defect = clean(C - Phi - jacobian(V,x)*f,1e-9);

    constraints = [
        sos(defect), ...
        C >= 0
        ];

    sol = solvesos(constraints,C,ops,[C; cV]);

    if sol.problem ~= 0
        warning('SOS status: %s', sol.info);
    end

    Cval = value(C);
end

function f = ks_galerkin_sdp(x,L)

    N = length(x);
    f = sdpvar(N,1);

    for k = 1:N

        q = k/L;

        lin = (q^2 - q^4)*x(k);

        nonlin = 0;

        for m = 1:N
            for n = 1:N
                coeff = triple_coeff(k,m,n,L);
                if coeff ~= 0
                    nonlin = nonlin + coeff*x(m)*x(n);
                end
            end
        end

        f(k) = lin + nonlin;
    end
end

function f = ks_galerkin_numeric(b,L)

    N = length(b);
    f = zeros(N,1);

    for k = 1:N

        q = k/L;
        f(k) = (q^2 - q^4)*b(k);

        for m = 1:N
            for n = 1:N
                coeff = triple_coeff(k,m,n,L);
                if coeff ~= 0
                    f(k) = f(k) + coeff*b(m)*b(n);
                end
            end
        end
    end
end

function c = triple_coeff(k,m,n,L)

    s = 0;

    if m+n == k
        s = s + 1;
    end

    if m-n == k
        s = s + 1;
    end

    if n-m == k
        s = s - 1;
    end

    c = -(n/(2*L))*s;
end

function exps = invariant_exponents_half_shift(N,dmax)

    exps = [];

    for d = 0:dmax

        E = exponent_rows(N,d);

        weights = E*(1:N).';
        keep = mod(weights,2) == 0;

        exps = [exps; E(keep,:)]; 
    end
end

function E = exponent_rows(n,d)

    if n == 1
        E = d;
    else
        E = [];
        for k = 0:d
            Ek = exponent_rows(n-1,d-k);
            E = [E; [k*ones(size(Ek,1),1), Ek]]; 
        end
    end
end

function mon = monomials_from_exponents(x,exps)

    nMon = size(exps,1);
    mon = sdpvar(nMon,1);

    for i = 1:nMon

        p = 1;

        for j = 1:length(x)
            if exps(i,j) > 0
                p = p*x(j)^exps(i,j);
            end
        end

        mon(i) = p;
    end
end

function [Bss, Ess] = find_steady_states_KS(N,L)

    opts = optimoptions('fsolve', ...
        'Display','off', ...
        'FunctionTolerance',1e-11, ...
        'StepTolerance',1e-11, ...
        'MaxIterations',1000, ...
        'MaxFunctionEvaluations',2e5);

    guesses = zeros(1,N);

    for k = 1:N

        q = k/L;

        if q^2 - q^4 > 0
            for amp = [-4 -2 -1 1 2 4]
                g = zeros(1,N);
                g(k) = amp;
                guesses = [guesses; g]; 
            end
        end
    end

    rng(1);

    for r = 1:40
        g = zeros(1,N);
        g(1:min(6,N)) = 3*randn(1,min(6,N));
        guesses = [guesses; g]; 
    end

    roots = [];

    for i = 1:size(guesses,1)

        g = guesses(i,:).';
        [b,~,exitflag] = fsolve(@(z) ks_galerkin_numeric(z,L),g,opts);

        if exitflag > 0 && norm(ks_galerkin_numeric(b,L)) < 1e-7

            if isempty(roots)
                roots = b.';
            else
                d = min(vecnorm(roots - b.',2,2));
                if d > 1e-5
                    roots = [roots; b.']; 
                end
            end
        end
    end

    Bss = roots;
    Ess = 0.5*sum(Bss.^2,2);

    [Ess,idx] = sort(Ess);
    Bss = Bss(idx,:);
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