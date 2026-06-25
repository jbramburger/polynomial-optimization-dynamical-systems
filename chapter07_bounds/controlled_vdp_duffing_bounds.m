%% controlled_vdp_duffing_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Controlled time-average bounds for coupled Van der Pol--Duffing oscillators
%
% This script synthesizes a polynomial feedback controller that reduces an
% upper bound on the long-time averaged synchronization error
%
%     Phi = (x1-x2)^2 + (y1-y2)^2
%
% for two mismatched coupled Van der Pol--Duffing oscillators.
%
% The SOS computations are performed in scaled variables z_i = x_i/R_i so
% that the computational domain is the unit box |z_i| <= 1.
%
% Pressing Run should reproduce the numerical output and save:
%
%     controlled_vdp_duffing_results.csv
%     controlled_vdp_duffing_t_x.pdf
%     controlled_vdp_duffing_x1_x2.pdf
%     controlled_vdp_duffing_bounds_by_degree.pdf
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

mu1   = 1;
mu2   = 1.2;
kappa = 0.01;

feedbackDeg  = 1;
relaxDegVals = [4 6 8];

maxAltIter = 20;
tolAlt     = 1e-4;

kCoeffMax = 10;

solver_name = 'mosek';
verbose_solver = 0;

saveFigures = true;
saveTable   = true;

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.black  = [0, 0, 0];
S.gray   = [0.40, 0.40, 0.40];

%% Build computational box from uncontrolled simulations

rng(1);

Tbox = 300;
tspanBox = linspace(0,Tbox,5000);
odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);

nBoxTraj = 20;
X0box = 3*(2*rand(4,nBoxTraj)-1);

fprintf('\n============================================================\n');
fprintf('Controlled Van der Pol--Duffing time-average bounds\n');
fprintf('============================================================\n\n');

fprintf('Simulating uncontrolled system to build computational box...\n');

Xall = [];

for j = 1:nBoxTraj
    [~,Xtmp] = ode45(@(t,x)vdp_duffing_rhs(t,x,mu1,mu2,kappa,0), ...
        tspanBox, X0box(:,j), odeopts);
    Xall = [Xall; Xtmp]; 
end

idxKeep = round(0.4*size(Xall,1)):size(Xall,1);
Xtail = Xall(idxKeep,:);

margin = 1.25;
R = margin*max(abs(Xtail),[],1).';
R = max(R,0.5);

fprintf('Physical box radii R = [%.3f %.3f %.3f %.3f]\n',R);
fprintf('SOS domain: |z_i| <= 1, with z_i = x_i/R_i.\n\n');

%% Initial feedback law in scaled variables

k0.coeffs = zeros(numOddMonomials(4,feedbackDeg),1);
k0.exps   = parityExponents(4,feedbackDeg,1);

%% Alternating SOS synthesis

Results = table();

bestC = Inf;
bestController = k0;
foundController = false;

kWarm = k0;
haveWarmStart = false;

for d = relaxDegVals

    fprintf('\n============================================================\n');
    fprintf('Relaxation degree d = %d\n',d);
    fprintf('============================================================\n');

    dk = feedbackDeg;

    dV = d + 1 - max(3,dk);

    if dV < 0
        warning('Relaxation degree too small for chosen feedback degree.');
        continue;
    end

    fprintf('Auxiliary degree dV = %d\n',dV);
    fprintf('Feedback degree dk  = %d\n',dk);

    if haveWarmStart
        fprintf('Using warm-start feedback from previous degree.\n');
    else
        fprintf('Using zero-feedback initialization.\n');
    end

    kCurrent = kWarm;
    Cold = Inf;

    bestDegreeC = Inf;
    bestDegreeController = kCurrent;
    bestDegreeIter = NaN;
    degreeSucceeded = false;
    usedWarmStartThisDegree = haveWarmStart;

    for iter = 1:maxAltIter

        fprintf('\nDegree %d, alternating iteration %d\n',d,iter);

        %% Auxiliary-function step

        [Vsol,Caux,solAux] = solve_aux_step_scaled( ...
            d,dV,dk,kCurrent,R,mu1,mu2,kappa,solver_name,verbose_solver);

        if solAux.problem ~= 0 || isnan(Caux) || isempty(Vsol.coeffs)
            warning('Auxiliary step failed: %s',solAux.info);
            break;
        end

        fprintf('  Auxiliary step bound C = %.8g\n',Caux);

        %% Feedback step

        [kNew,Cnew,solK] = solve_feedback_step_scaled( ...
            d,dV,dk,Vsol,R,mu1,mu2,kappa,kCoeffMax,solver_name,verbose_solver);

        if solK.problem ~= 0 || isnan(Cnew) || isempty(kNew.coeffs)
            warning('Feedback step failed: %s',solK.info);
            break;
        end

        kInf = norm(kNew.coeffs,inf);
        kTwo = norm(kNew.coeffs,2);

        fprintf('  Feedback step bound  C = %.8g\n',Cnew);
        fprintf('  ||k||_inf = %.3e, ||k||_2 = %.3e\n', kInf, kTwo);

        Results = [Results; table(d,iter,Caux,Cnew,kInf,kTwo, ...
            usedWarmStartThisDegree,solAux.problem,solK.problem, ...
            'VariableNames',{'degree','iteration','C_aux','C_feedback', ...
            'k_inf','k_2','used_warm_start','status_aux','status_feedback'})]; 

        if Cnew < bestDegreeC
            bestDegreeC = Cnew;
            bestDegreeController = kNew;
            bestDegreeIter = iter;
            degreeSucceeded = true;
        end

        if Cnew < bestC
            bestC = Cnew;
            bestController = kNew;
            foundController = true;
        end

        if isfinite(Cold) && abs(Cnew-Cold) < tolAlt
            fprintf('  Alternation converged: |Delta C| = %.3e\n',abs(Cnew-Cold));
            break;
        end

        Cold = Cnew;
        kCurrent = kNew;
    end

    if degreeSucceeded
        kWarm = bestDegreeController;
        haveWarmStart = true;

        fprintf('\nBest controller at degree %d came from iteration %d with C = %.8g.\n', ...
            d,bestDegreeIter,bestDegreeC);
    else
        fprintf('\nNo successful controller found at degree %d.\n',d);
    end
end

disp(Results);

if saveTable && ~isempty(Results)
    writetable(Results,'controlled_vdp_duffing_results.csv');
    fprintf('Saved table: controlled_vdp_duffing_results.csv\n');
end

if foundController
    fprintf('\nBest certified upper bound C = %.10g\n',bestC);
    fprintf('Best feedback coefficients in scaled variables z=x./R:\n');
    disp(bestController.coeffs.');
else
    warning('No successful feedback step completed. Using zero feedback for plots.');
    bestC = NaN;
    bestController = k0;
end

%% Compare uncontrolled and controlled simulations

Tplot = 120;
tspan = linspace(0,Tplot,3000);

x0 = [2.0; 0.2; -1.5; 0.1];

[~,Xu] = ode45(@(t,x)vdp_duffing_rhs(t,x,mu1,mu2,kappa,0), ...
    tspan, x0, odeopts);

[~,Xc] = ode45(@(t,x)vdp_duffing_rhs(t,x,mu1,mu2,kappa, ...
    eval_feedback_physical(bestController,x,R)), ...
    tspan, x0, odeopts);

Phi_u = (Xu(:,1)-Xu(:,3)).^2 + (Xu(:,2)-Xu(:,4)).^2;
Phi_c = (Xc(:,1)-Xc(:,3)).^2 + (Xc(:,2)-Xc(:,4)).^2;

avgPhi_u = cumtrapz(tspan,Phi_u)./max(tspan(:),eps);
avgPhi_c = cumtrapz(tspan,Phi_c)./max(tspan(:),eps);

fprintf('\nEmpirical uncontrolled average at final time = %.8g\n',avgPhi_u(end));
fprintf('Empirical controlled average at final time   = %.8g\n',avgPhi_c(end));

%% Figure 1: time traces

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 13 8]);

hold on
box on

plot(tspan,Xu(:,1), '--', ...
    'Color',S.gray, ...
    'LineWidth',2.0);

plot(tspan,Xc(:,3), ...
    'Color',S.red, ...
    'LineWidth',4);

plot(tspan,Xc(:,1), ...
    'Color',S.blue, ...
    'LineWidth',4);

xlabel('$t$','Interpreter','latex','FontSize',36);

legend({'uncontrolled $x_1$', 'controlled $x_2$', 'controlled $x_1$'}, ...
    'Interpreter','latex','Location','northeast','FontSize',20);

set(gca,'FontSize',24,'TickLabelInterpreter','latex','LineWidth',1.0);
xlim([0 40])
grid on

if saveFigures
    save_pdf_figure(fig1,'controlled_vdp_duffing_t_x.pdf');
end

%% Figure 2: synchronization projection

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 11 11]);

hold on
box on

plot(Xu(:,1),Xu(:,3), ...
    'Color',S.gray, ...
    'LineWidth',1);

rr = max(R([1 3]));
plot([-rr rr],[-rr rr], ...
    '--','Color',S.black,'LineWidth',5);

plot(Xc(:,1),Xc(:,3), ...
    'Color',S.blue, ...
    'LineWidth',5);

xlabel('$x_1$','Interpreter','latex','FontSize',36);
ylabel('$x_2$','Interpreter','latex','FontSize',36);

legend({'uncontrolled','$x_1=x_2$','controlled'}, ...
    'Interpreter','latex','Location','best','FontSize',24);

set(gca,'FontSize',24,'TickLabelInterpreter','latex','LineWidth',1.0);
xlim([-2 2])
ylim([-2 2])
grid on

if saveFigures
    save_pdf_figure(fig2,'controlled_vdp_duffing_x1_x2.pdf');
end

%% Figure 3: certified bounds by degree

if ~isempty(Results)

    fig3 = figure;
    set(fig3,'Color','w','Units','centimeters','Position',[2 2 11 8]);

    hold on
    box on

    degList = unique(Results.degree);
    bestByDeg = nan(size(degList));

    for i = 1:length(degList)
        rows = Results.degree == degList(i);
        bestByDeg(i) = min(Results.C_feedback(rows),[],'omitnan');
    end

    plot(degList,bestByDeg,'-o', ...
        'Color',S.blue, ...
        'LineWidth',2.8, ...
        'MarkerSize',8);

    xlabel('relaxation degree $d$','Interpreter','latex','FontSize',28);
    ylabel('best certified upper bound $C_d$', ...
        'Interpreter','latex','FontSize',28);

    set(gca,'FontSize',24,'TickLabelInterpreter','latex','LineWidth',1.0);
    grid on

    if saveFigures
        save_pdf_figure(fig3,'controlled_vdp_duffing_bounds_by_degree.pdf');
    end
end

fprintf('\nFinished controlled Van der Pol--Duffing example.\n');

%% Local functions

function [Vsol,Cval,sol] = solve_aux_step_scaled( ...
    d,dV,dk,kFixed,R,mu1,mu2,kappa,solver_name,verbose)

    yalmip('clear');

    z = sdpvar(4,1);
    C = sdpvar(1,1);

    xphys = R(:).*z;

    fphys0 = vdp_duffing_sdp(xphys,mu1,mu2,kappa,0);
    f0 = fphys0./R(:);

    g = [0; 1/R(2); 0; 0];

    kExpr = polyFromExps(z,kFixed.coeffs,kFixed.exps);

    Phi = (xphys(1)-xphys(3))^2 + (xphys(2)-xphys(4))^2;

    expsV = parityExponents(4,dV,0);
    cV = sdpvar(size(expsV,1),1);

    idxConst = find(sum(expsV,2) == 0);

    consGauge = [];
    if ~isempty(idxConst)
        consGauge = [cV(idxConst) == 0];
    end

    zV = monomialVectorFromExps(z,expsV);
    V = cV.'*zV;

    D = clean(C - Phi - jacobian(V,z)*(f0 + g*kExpr),1e-10);

    [cons,params] = sos_unit_box_certificate(D,z,d);
    cons = [cons, consGauge];

    params = [C; cV; params];

    opts = solver_opts(solver_name,verbose);

    sol = solvesos(cons,C,opts,params);

    Vsol.exps = expsV;

    if sol.problem == 0

        Cval = value(C);
        Vcoeffs = value(cV);

        if ~isempty(idxConst) && isnan(Vcoeffs(idxConst))
            Vcoeffs(idxConst) = 0;
        end

        if isnan(Cval) || any(isnan(Vcoeffs))
            sol.problem = 99;
            sol.info = 'NaN returned in auxiliary solution';
            Cval = NaN;
            Vsol.coeffs = [];
        else
            Vsol.coeffs = Vcoeffs;
        end

    else

        Cval = NaN;
        Vsol.coeffs = [];
    end
end

function [ksol,Cval,sol] = solve_feedback_step_scaled( ...
    d,dV,dk,Vfixed,R,mu1,mu2,kappa,kCoeffMax,solver_name,verbose)

    yalmip('clear');

    z = sdpvar(4,1);
    C = sdpvar(1,1);

    xphys = R(:).*z;

    fphys0 = vdp_duffing_sdp(xphys,mu1,mu2,kappa,0);
    f0 = fphys0./R(:);

    g = [0; 1/R(2); 0; 0];

    expsK = parityExponents(4,dk,1);
    cK = sdpvar(size(expsK,1),1);
    zK = monomialVectorFromExps(z,expsK);
    kExpr = cK.'*zK;

    V = polyFromExps(z,Vfixed.coeffs,Vfixed.exps);

    Phi = (xphys(1)-xphys(3))^2 + (xphys(2)-xphys(4))^2;

    D = clean(C - Phi - jacobian(V,z)*(f0 + g*kExpr),1e-10);

    [cons,params] = sos_unit_box_certificate(D,z,d);

    cons = [cons, -kCoeffMax <= cK <= kCoeffMax];

    params = [C; cK; params];

    opts = solver_opts(solver_name,verbose);

    sol = solvesos(cons,C,opts,params);

    ksol.exps = expsK;

    if sol.problem == 0

        Cval = value(C);
        Kcoeffs = value(cK);

        if isnan(Cval) || any(isnan(Kcoeffs))
            sol.problem = 99;
            sol.info = 'NaN returned in feedback solution';
            Cval = NaN;
            ksol.coeffs = [];
        else
            ksol.coeffs = Kcoeffs;
        end

    else

        Cval = NaN;
        ksol.coeffs = [];
    end
end

function [cons,params] = sos_unit_box_certificate(D,z,d)

    n = length(z);
    cons = [];
    params = [];

    [s0,c0] = polynomial(z,d);
    cons = [cons, sos(s0)];
    params = [params; c0];

    rhs = s0;

    multDeg = max(d - 2,0);

    for j = 1:n

        gj = 1 - z(j)^2;

        [sj,cj] = polynomial(z,multDeg);

        cons = [cons, sos(sj)]; 
        params = [params; cj]; 

        rhs = rhs + sj*gj;
    end

    p = clean(D - rhs,1e-10);
    coef = coefficients(p,z);

    cons = [cons, coef(:) == 0];
end

function opts = solver_opts(solver_name,verbose)

    opts = sdpsettings('solver',solver_name,'verbose',verbose);

    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-7;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-7;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-7;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-7;
end

function dx = vdp_duffing_rhs(~,x,mu1,mu2,kappa,u)

    dx = zeros(4,1);

    dx(1) = x(2);

    dx(2) = x(1) - x(1)^3 ...
        + mu1*(1 - x(1)^2)*x(2) ...
        + kappa*(x(3)-x(1)) + u;

    dx(3) = x(4);

    dx(4) = x(3) - x(3)^3 ...
        + mu2*(1 - x(3)^2)*x(4) ...
        + kappa*(x(1)-x(3));
end

function f = vdp_duffing_sdp(x,mu1,mu2,kappa,u)

    f = sdpvar(4,1);

    f(1) = x(2);

    f(2) = x(1) - x(1)^3 ...
        + mu1*(1 - x(1)^2)*x(2) ...
        + kappa*(x(3)-x(1)) + u;

    f(3) = x(4);

    f(4) = x(3) - x(3)^3 ...
        + mu2*(1 - x(3)^2)*x(4) ...
        + kappa*(x(1)-x(3));
end

function u = eval_feedback_physical(controller,xphys,R)

    z = xphys(:)./R(:);
    u = evalPoly(controller.coeffs,controller.exps,z);
end

function nmon = numOddMonomials(n,maxDeg)

    exps = parityExponents(n,maxDeg,1);
    nmon = size(exps,1);
end

function exps = parityExponents(n,maxDeg,parity)

    allExps = generateExponents(n,maxDeg);
    totalDeg = sum(allExps,2);

    exps = allExps(mod(totalDeg,2)==parity,:);
end

function z = monomialVectorFromExps(x,exps)

    z = [];

    for k = 1:size(exps,1)

        mon = 1;

        for j = 1:length(x)
            if exps(k,j) > 0
                mon = mon*x(j)^exps(k,j);
            end
        end

        z = [z; mon]; 
    end
end

function p = polyFromExps(x,coeffs,exps)

    if isempty(coeffs) || any(isnan(coeffs))
        error('polyFromExps:InvalidCoefficients', ...
            'Polynomial coefficients are empty or contain NaNs.');
    end

    z = monomialVectorFromExps(x,exps);
    p = coeffs(:).'*z;
end

function val = evalPoly(coeffs,exps,x)

    x = x(:);
    val = 0;

    if isempty(coeffs)
        return;
    end

    for k = 1:size(exps,1)

        mon = 1;

        for j = 1:length(x)
            mon = mon*x(j)^exps(k,j);
        end

        val = val + coeffs(k)*mon;
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