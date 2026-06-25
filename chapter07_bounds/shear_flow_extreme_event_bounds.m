%% shear_flow_extreme_event_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Finite-time energy amplification in a shear-flow model
%
% This script computes finite-time auxiliary-function bounds for transient
% energy growth in the shifted MFE9 shear-flow model.
%
% Pressing Run should reproduce the numerical output and save:
%
%     shear_extreme_results.csv
%     shear_extreme.pdf
%     shear_extreme_bounds.pdf
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

Re = 20;
Lx = 4*pi;
Lz = 2*pi;

T  = 20;
E0 = 0.1;
EK = 0.5;

spaceDegVals = [2 4 6];

timeDegVals_quad    = [2 4 6 8];
timeDegVals_quartic = [2 4];
timeDegVals_sextic  = [];

NsGrid = 61;

% Use 100000 to reproduce the high-resolution Monte Carlo search used for
% the book figure. The default is smaller so the script runs quickly.
nTraj = 2000;

verbose_solver = 0;

saveFigures = true;
saveTable = true;

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.black  = [0, 0, 0];
S.gray   = [0.40, 0.40, 0.40];
S.light_gray = [0.78, 0.78, 0.78];

%% Monte Carlo simulations

rng(1);

tspan = linspace(0,T,1200);
odeopts = odeset('RelTol',1e-9,'AbsTol',1e-11);

X0 = zeros(9,nTraj);

for j = 1:nTraj
    v = randn(9,1);
    v = v/norm(v);
    X0(:,j) = sqrt(2*E0)*v;
end

Etraj = nan(length(tspan),nTraj);
maxE = nan(nTraj,1);

fprintf('\n============================================================\n');
fprintf('Finite-time energy amplification in the shifted MFE9 model\n');
fprintf('============================================================\n\n');

fprintf('Re = %.6g\n', Re);
fprintf('T  = %.6g\n', T);
fprintf('E0 = %.6g\n', E0);
fprintf('EK = %.6g\n', EK);
fprintf('Monte Carlo trajectories: %d\n\n', nTraj);

fprintf('Simulating sample trajectories...\n');

for j = 1:nTraj

    [~,xout] = ode45(@(t,x)mfe9_rhs_shifted(t,x,Re,Lx,Lz), ...
        tspan, X0(:,j), odeopts);

    E = 0.5*sum(xout.^2,2);

    Etraj(:,j) = E;
    maxE(j) = max(E);
end

[bestSampleE,idxBest] = max(maxE);

fprintf('\nLargest sampled energy = %.8g\n', bestSampleE);
fprintf('Computational energy ball EK = %.8g\n', EK);

%% Solve auxiliary-function bounds

allResults = [];

for sd = spaceDegVals

    if sd == 2
        timeDegVals = timeDegVals_quad;
    elseif sd == 4
        timeDegVals = timeDegVals_quartic;
    else
        timeDegVals = timeDegVals_sextic;
    end

    for td = timeDegVals

        fprintf('\nSolving space degree %d, time degree %d...\n', sd, td);

        if sd == 2
            [Ctmp, sol] = solve_mfe9_extreme_sdp_sym_quadV( ...
                Re,Lx,Lz,T,E0,EK,td,NsGrid,verbose_solver);
        else
            [Ctmp, sol] = solve_mfe9_extreme_sos_sym_polyV( ...
                Re,Lx,Lz,T,E0,EK,td,sd,verbose_solver);
        end

        if sol.problem == 0
            fprintf('  Bound = %.8g\n', Ctmp);
        else
            Ctmp = NaN;
            warning('Solver issue: %s', sol.info);
        end

        allResults = [allResults; sd, td, Ctmp, sol.problem]; 
    end
end

Results = array2table(allResults, ...
    'VariableNames', {'space_degree','time_degree','upper_bound','solver_status'});

disp(Results);

bestBound = min(Results.upper_bound,[],'omitnan');

fprintf('\nLargest sampled energy = %.8g\n', bestSampleE);
fprintf('Best computed bound    = %.8g\n', bestBound);

if saveTable
    writetable(Results,'shear_extreme_results.csv');
    fprintf('Saved table: shear_extreme_results.csv\n');
end

%% Figure 1: energy-time plot

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 13 8]);

hold on
box on

plotStride = max(1,floor(nTraj/500));

for j = 1:plotStride:nTraj
    plot(tspan,Etraj(:,j), ...
        'Color',0.82*S.light_gray, ...
        'LineWidth',2);
end

plot(tspan,Etraj(:,idxBest), ...
    'Color',S.red, ...
    'LineWidth',5);

yline(E0,'--','Color',S.blue,'LineWidth',5.0);

if ~isnan(bestBound)
    yline(bestBound,'-','Color',S.black,'LineWidth',5);
end

xlabel('$t$','Interpreter','latex','FontSize',36,'FontWeight','bold');
ylabel('$E(t)=\frac12\|x(t)\|^2$', ...
    'Interpreter','latex','FontSize',36,'FontWeight','bold');

set(gca,'FontSize',36,'TickLabelInterpreter','latex','LineWidth',1.0);
grid on

if saveFigures
    save_pdf_figure(fig1,'shear_extreme.pdf');
end

%% Figure 2: bound comparison

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 12 9]);

hold on
box on

for sd = spaceDegVals

    rows = Results.space_degree == sd & ~isnan(Results.upper_bound);

    if any(rows)
        plot(Results.time_degree(rows),Results.upper_bound(rows),'-o', ...
            'LineWidth',2.4, ...
            'MarkerSize',8, ...
            'DisplayName',sprintf('space degree %d',sd));
    end
end

yline(bestSampleE,'--','Color',S.red,'LineWidth',2.0, ...
    'DisplayName','largest sampled value');

xlabel('time polynomial degree','Interpreter','latex','FontSize',22);
ylabel('energy bound','Interpreter','latex','FontSize',22);

legend('Interpreter','latex','Location','best','FontSize',18);

set(gca,'FontSize',20,'TickLabelInterpreter','latex','LineWidth',1.0);
grid on

if saveFigures
    save_pdf_figure(fig2,'shear_extreme_bounds.pdf');
end

fprintf('\nFinished shear-flow extreme-event example.\n');

%% Local functions

function [Cval, sol] = solve_mfe9_extreme_sdp_sym_quadV( ...
    Re,Lx,Lz,T,E0,EK,timeDeg,NsGrid,verbose)

    yalmip('clear');

    n = 9;
    x = sdpvar(n,1);
    C = sdpvar(1,1);

    f = mfe9_rhs_shifted_sdp(x,Re,Lx,Lz);

    b = replace(f,x,zeros(n,1));
    J = replace(jacobian(f,x),x,zeros(n,1));
    N = clean(f - b - J*x,1e-12);

    A = 0.5*(double(J) + double(J).');

    cancellationResidual = test_energy_cancellation(A,Re,Lx,Lz);
    fprintf('  nonlinear energy cancellation residual = %.3e\n', cancellationResidual);

    G = [
         1 -1 -1 -1 -1  1  1  1  1
         1  1  1 -1 -1 -1 -1 -1  1
        ];

    Pcell = cell(timeDeg+1,1);
    rcell = cell(timeDeg+1,1);
    gcell = sdpvar(timeDeg+1,1);

    cons = [];

    for ell = 0:timeDeg

        Pell = sdpvar(n,n,'symmetric');
        rell = sdpvar(n,1);

        for i = 1:n

            for j = i:n
                allowed = all(G(:,i).*G(:,j) == 1);
                if ~allowed
                    cons = [cons, Pell(i,j) == 0]; 
                end
            end

            allowedLinear = all(G(:,i) == 1);

            if ~allowedLinear
                cons = [cons, rell(i) == 0]; 
            end
        end

        Pcell{ell+1} = Pell;
        rcell{ell+1} = rell;

        cubic = clean(x.'*Pell*N,1e-10);
        cc = coefficients(cubic(:),x);

        cons = [cons, cc(:) == 0]; 
    end

    sgrid = linspace(0,1,NsGrid);

    for q = 1:NsGrid

        ss = sgrid(q);

        [P,Ps,r,rs,gam,gams] = eval_aux_coeffs( ...
            ss,Pcell,rcell,gcell,timeDeg);

        V  = 0.5*x.'*P*x + r.'*x + gam;
        Vs = 0.5*x.'*Ps*x + rs.'*x + gams;

        LieV = jacobian(0.5*x.'*P*x + r.'*x,x)*(T*(b + J*x)) ...
             + r.'*(T*N);

        D1 = clean(-Vs - LieV,1e-10);
        D2 = clean(V - 0.5*(x.'*x),1e-10);

        cons = add_ball_sprocedure(cons,D1,x,EK);
        cons = add_global_quadratic_psd(cons,D2,x);
    end

    P0 = Pcell{1};
    r0 = rcell{1};
    g0 = gcell(1);

    V0 = 0.5*x.'*P0*x + r0.'*x + g0;
    D3 = clean(C - V0,1e-10);

    cons = add_ball_sprocedure(cons,D3,x,E0);
    cons = [cons, C >= 0];

    opts = sdpsettings('solver','mosek','verbose',verbose);

    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-8;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-8;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-8;

    sol = optimize(cons,C,opts);

    if sol.problem == 0
        Cval = value(C);
    else
        Cval = NaN;
    end
end

function [Cval, sol] = solve_mfe9_extreme_sos_sym_polyV( ...
    Re,Lx,Lz,T,E0,EK,timeDeg,spaceDeg,verbose)

    yalmip('clear');

    n = 9;
    x = sdpvar(n,1);
    s = sdpvar(1,1);
    C = sdpvar(1,1);

    vars_xs = [x; s];

    f = mfe9_rhs_shifted_sdp(x,Re,Lx,Lz);

    E = 0.5*(x.'*x);

    gK = EK - E;
    gX = E0 - E;
    gS = s*(1-s);

    G = [
         1 -1 -1 -1 -1  1  1  1  1
         1  1  1 -1 -1 -1 -1 -1  1
        ];

    z = invariant_monomial_basis(x,spaceDeg,G);

    V = 0;
    coeffsV = [];

    for ell = 0:timeDeg

        c = sdpvar(length(z),1);

        V = V + s^ell*(c.'*z);
        coeffsV = [coeffsV; c]; 
    end

    Vdot = jacobian(V,x)*(T*f) + jacobian(V,s);

    D1 = clean(-Vdot,1e-10);
    D2 = clean(V - E,1e-10);
    D3 = clean(C - replace(V,s,0),1e-10);

    if spaceDeg == 4
        multDegK = 4;
        multDegS = 4;
        multDegX = 2;
    elseif spaceDeg == 6
        multDegK = 6;
        multDegS = 6;
        multDegX = 4;
    else
        multDegK = spaceDeg;
        multDegS = spaceDeg;
        multDegX = max(spaceDeg-2,0);
    end

    [m1K,c1K] = polynomial(vars_xs,multDegK);
    [m1S,c1S] = polynomial(vars_xs,multDegS);
    [m3X,c3X] = polynomial(x,multDegX);

    cons = [
        sos(m1K), ...
        sos(m1S), ...
        sos(m3X), ...
        sos(D1 - m1K*gK - m1S*gS), ...
        sos(D2), ...
        sos(D3 - m3X*gX), ...
        C >= 0
        ];

    params = [C; coeffsV; c1K; c1S; c3X];

    opts = sdpsettings('solver','mosek','verbose',verbose);

    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-7;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-7;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-7;

    sol = solvesos(cons,C,opts,params);

    if sol.problem == 0
        Cval = value(C);
    else
        Cval = NaN;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Monomial basis helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = invariantMonomialBasis(x,maxDeg,G)

    n = length(x);
    exps = generateExponents(n,maxDeg);

    keep = false(size(exps,1),1);

    for k = 1:size(exps,1)
        alpha = exps(k,:);

        invariant = true;
        for g = 1:size(G,1)
            parity = mod(sum(alpha .* (G(g,:) == -1)),2);
            if parity ~= 0
                invariant = false;
                break;
            end
        end

        keep(k) = invariant;
    end

    exps = exps(keep,:);

    z = [];
    for k = 1:size(exps,1)
        mon = 1;
        for j = 1:n
            if exps(k,j) > 0
                mon = mon*x(j)^exps(k,j);
            end
        end
        z = [z; mon];
    end
end

function exps = generateExponents(n,maxDeg)

    exps = [];
    current = zeros(1,n);

    for totalDeg = 0:maxDeg
        exps = [exps; genFixedDegree(n,totalDeg,current,1)];
    end
end

function exps = genFixedDegree(n,totalDeg,current,idx)

    if idx == n
        current(idx) = totalDeg;
        exps = current;
        return;
    end

    exps = [];

    for a = 0:totalDeg
        current(idx) = a;
        exps = [exps; genFixedDegree(n,totalDeg-a,current,idx+1)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quadratic auxiliary-function helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,Ps,r,rs,gam,gams] = eval_aux_coeffs( ...
    s,Pcell,rcell,gcell,timeDeg)

    P  = 0*Pcell{1};
    Ps = 0*Pcell{1};

    r  = 0*rcell{1};
    rs = 0*rcell{1};

    gam  = 0*gcell(1);
    gams = 0*gcell(1);

    for ell = 0:timeDeg

        P   = P   + s^ell*Pcell{ell+1};
        r   = r   + s^ell*rcell{ell+1};
        gam = gam + s^ell*gcell(ell+1);

        if ell >= 1
            Ps   = Ps   + ell*s^(ell-1)*Pcell{ell+1};
            rs   = rs   + ell*s^(ell-1)*rcell{ell+1};
            gams = gams + ell*s^(ell-1)*gcell(ell+1);
        end
    end
end

function cons = addBallSprocedure(cons,q,x,Emax)

    n = length(x);
    mu = sdpvar(1,1);

    qcert = clean(q - mu*(Emax - 0.5*(x.'*x)),1e-10);

    H = jacobian(jacobian(qcert,x).',x);
    Q = 0.5*H;

    l = replace(jacobian(qcert,x).',x,zeros(n,1));
    c = replace(qcert,x,zeros(n,1));

    M = [c, 0.5*l.'; 0.5*l, Q];

    cons = [cons, mu >= 0, M >= 0];
end

function cons = addGlobalQuadraticPSD(cons,q,x)

    n = length(x);

    q = clean(q,1e-10);

    H = jacobian(jacobian(q,x).',x);
    Q = 0.5*H;

    l = replace(jacobian(q,x).',x,zeros(n,1));
    c = replace(q,x,zeros(n,1));

    M = [c, 0.5*l.'; 0.5*l, Q];

    cons = [cons, M >= 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residual = test_energy_cancellation(A,Re,Lx,Lz)

    rng(10);
    residual = 0;

    for k = 1:200
        x = randn(9,1);
        x = x/norm(x);

        f = mfe9_rhs_shifted(0,x,Re,Lx,Lz);

        q_true = x.'*f;
        q_quad = x.'*A*x;

        residual = max(residual,abs(q_true-q_quad));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shifted MFE9 model: symbolic and numeric versions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = mfe9_rhs_shifted_sdp(x,Re,Lx,Lz)

    a = [x(1)+1;
         x(2);
         x(3);
         x(4);
         x(5);
         x(6);
         x(7);
         x(8);
         x(9)];

    f = mfe9_rhs_a_sdp(a,Re,Lx,Lz);
end

function dx = mfe9_rhs_shifted(~,x,Re,Lx,Lz)

    a = [x(1)+1;
         x(2);
         x(3);
         x(4);
         x(5);
         x(6);
         x(7);
         x(8);
         x(9)];

    dx = mfe9_rhs_a(a,Re,Lx,Lz);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFE9 equations in original a-coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function da = mfe9_rhs_a_sdp(a,Re,Lx,Lz)

    alpha = 2*pi/Lx;
    beta  = pi/2;
    gamma = 2*pi/Lz;

    kag  = sqrt(alpha^2 + gamma^2);
    kbg  = sqrt(beta^2  + gamma^2);
    kabg = sqrt(alpha^2 + beta^2 + gamma^2);

    da = sdpvar(9,1);

    da(1) = beta^2/Re*(1-a(1)) ...
        - sqrt(3/2)*beta*gamma/kabg*a(6)*a(8) ...
        + sqrt(3/2)*beta*gamma/kbg*a(2)*a(3);

    da(2) = -((4/3)*beta^2 + gamma^2)/Re*a(2) ...
        + (5*sqrt(2)/(3*sqrt(3)))*gamma^2/kag*a(4)*a(6) ...
        - gamma^2/(sqrt(6)*kag)*a(5)*a(7) ...
        - alpha*beta*gamma/(sqrt(6)*kag*kabg)*a(5)*a(8) ...
        - sqrt(3/2)*beta*gamma/kbg*a(1)*a(3) ...
        - sqrt(3/2)*beta*gamma/kbg*a(3)*a(9);

    da(3) = -(beta^2 + gamma^2)/Re*a(3) ...
        + (2/sqrt(6))*alpha*beta*gamma/(kag*kbg)*(a(4)*a(7)+a(5)*a(6)) ...
        + (beta^2*(3*alpha^2+gamma^2)-3*gamma^2*(alpha^2+gamma^2)) ...
          /(sqrt(6)*kag*kbg*kabg)*a(4)*a(8);

    da(4) = -(3*alpha^2 + 4*beta^2)/(3*Re)*a(4) ...
        - alpha/sqrt(6)*a(1)*a(5) ...
        - (10/(3*sqrt(6)))*alpha^2/kag*a(2)*a(6) ...
        - sqrt(3/2)*alpha*beta*gamma/(kag*kbg)*a(3)*a(7) ...
        - sqrt(3/2)*alpha^2*beta^2/(kag*kbg*kabg)*a(3)*a(8) ...
        - alpha/sqrt(6)*a(5)*a(9);

    da(5) = -(alpha^2 + beta^2)/Re*a(5) ...
        + alpha/sqrt(6)*a(1)*a(4) ...
        + alpha^2/(sqrt(6)*kag)*a(2)*a(7) ...
        - alpha*beta*gamma/(sqrt(6)*kag*kabg)*a(2)*a(8) ...
        + alpha/sqrt(6)*a(4)*a(9) ...
        + (2/sqrt(6))*alpha*beta*gamma/(kag*kbg)*a(3)*a(6);

    da(6) = -(3*alpha^2 + 4*beta^2 + 3*gamma^2)/(3*Re)*a(6) ...
        + alpha/sqrt(6)*a(1)*a(7) ...
        + sqrt(3/2)*beta*gamma/kabg*a(1)*a(8) ...
        + (10/(3*sqrt(6)))*(alpha^2-gamma^2)/kag*a(2)*a(4) ...
        - 2*sqrt(2/3)*alpha*beta*gamma/(kag*kbg)*a(3)*a(5) ...
        + alpha/sqrt(6)*a(7)*a(9) ...
        + sqrt(3/2)*beta*gamma/kabg*a(8)*a(9);

    da(7) = -(alpha^2 + beta^2 + gamma^2)/Re*a(7) ...
        - alpha/sqrt(6)*(a(1)*a(6)+a(6)*a(9)) ...
        + (gamma^2-alpha^2)/(sqrt(6)*kag)*a(2)*a(5) ...
        + alpha*beta*gamma/(sqrt(6)*kag*kbg)*a(3)*a(4);

    da(8) = -(alpha^2 + beta^2 + gamma^2)/Re*a(8) ...
        + (2/sqrt(6))*alpha*beta*gamma/(kag*kabg)*a(2)*a(5) ...
        + gamma^2*(3*alpha^2-beta^2+3*gamma^2) ...
          /(sqrt(6)*kag*kbg*kabg)*a(3)*a(4);

    da(9) = -9*beta^2/Re*a(9) ...
        + sqrt(3/2)*beta*gamma/kbg*a(2)*a(3) ...
        - sqrt(3/2)*beta*gamma/kabg*a(6)*a(8);
end

function da = mfe9_rhs_a(a,Re,Lx,Lz)

    alpha = 2*pi/Lx;
    beta  = pi/2;
    gamma = 2*pi/Lz;

    kag  = sqrt(alpha^2 + gamma^2);
    kbg  = sqrt(beta^2  + gamma^2);
    kabg = sqrt(alpha^2 + beta^2 + gamma^2);

    da = zeros(9,1);

    da(1) = beta^2/Re*(1-a(1)) ...
        - sqrt(3/2)*beta*gamma/kabg*a(6)*a(8) ...
        + sqrt(3/2)*beta*gamma/kbg*a(2)*a(3);

    da(2) = -((4/3)*beta^2 + gamma^2)/Re*a(2) ...
        + (5*sqrt(2)/(3*sqrt(3)))*gamma^2/kag*a(4)*a(6) ...
        - gamma^2/(sqrt(6)*kag)*a(5)*a(7) ...
        - alpha*beta*gamma/(sqrt(6)*kag*kabg)*a(5)*a(8) ...
        - sqrt(3/2)*beta*gamma/kbg*a(1)*a(3) ...
        - sqrt(3/2)*beta*gamma/kbg*a(3)*a(9);

    da(3) = -(beta^2 + gamma^2)/Re*a(3) ...
        + (2/sqrt(6))*alpha*beta*gamma/(kag*kbg)*(a(4)*a(7)+a(5)*a(6)) ...
        + (beta^2*(3*alpha^2+gamma^2)-3*gamma^2*(alpha^2+gamma^2)) ...
          /(sqrt(6)*kag*kbg*kabg)*a(4)*a(8);

    da(4) = -(3*alpha^2 + 4*beta^2)/(3*Re)*a(4) ...
        - alpha/sqrt(6)*a(1)*a(5) ...
        - (10/(3*sqrt(6)))*alpha^2/kag*a(2)*a(6) ...
        - sqrt(3/2)*alpha*beta*gamma/(kag*kbg)*a(3)*a(7) ...
        - sqrt(3/2)*alpha^2*beta^2/(kag*kbg*kabg)*a(3)*a(8) ...
        - alpha/sqrt(6)*a(5)*a(9);

    da(5) = -(alpha^2 + beta^2)/Re*a(5) ...
        + alpha/sqrt(6)*a(1)*a(4) ...
        + alpha^2/(sqrt(6)*kag)*a(2)*a(7) ...
        - alpha*beta*gamma/(sqrt(6)*kag*kabg)*a(2)*a(8) ...
        + alpha/sqrt(6)*a(4)*a(9) ...
        + (2/sqrt(6))*alpha*beta*gamma/(kag*kbg)*a(3)*a(6);

    da(6) = -(3*alpha^2 + 4*beta^2 + 3*gamma^2)/(3*Re)*a(6) ...
        + alpha/sqrt(6)*a(1)*a(7) ...
        + sqrt(3/2)*beta*gamma/kabg*a(1)*a(8) ...
        + (10/(3*sqrt(6)))*(alpha^2-gamma^2)/kag*a(2)*a(4) ...
        - 2*sqrt(2/3)*alpha*beta*gamma/(kag*kbg)*a(3)*a(5) ...
        + alpha/sqrt(6)*a(7)*a(9) ...
        + sqrt(3/2)*beta*gamma/kabg*a(8)*a(9);

    da(7) = -(alpha^2 + beta^2 + gamma^2)/Re*a(7) ...
        - alpha/sqrt(6)*(a(1)*a(6)+a(6)*a(9)) ...
        + (gamma^2-alpha^2)/(sqrt(6)*kag)*a(2)*a(5) ...
        + alpha*beta*gamma/(sqrt(6)*kag*kbg)*a(3)*a(4);

    da(8) = -(alpha^2 + beta^2 + gamma^2)/Re*a(8) ...
        + (2/sqrt(6))*alpha*beta*gamma/(kag*kabg)*a(2)*a(5) ...
        + gamma^2*(3*alpha^2-beta^2+3*gamma^2) ...
          /(sqrt(6)*kag*kbg*kabg)*a(3)*a(4);

    da(9) = -9*beta^2/Re*a(9) ...
        + sqrt(3/2)*beta*gamma/kbg*a(2)*a(3) ...
        - sqrt(3/2)*beta*gamma/kabg*a(6)*a(8);
end