%% hybrid_clutch_control.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Extensions and Frontiers
% Example: SOS control synthesis for a hybrid clutch model
%
% This script synthesizes polynomial feedback controllers for a hybrid dry
% clutch model with open, slipping, and locked modes.
%
% The SOS problem prescribes an exponential decay rate alpha0 and searches
% for minimum-gain linear feedback laws satisfying
%
%     dV/dt <= -alpha0 V
%
% in the slipping and locked modes. The script also checks compatibility of
% the reset map from slipping to locked dynamics using the same SOS
% framework.
%
% Pressing Run should reproduce the numerical output and save:
%
%     hybrid_control_1.pdf
%     hybrid_control_2.pdf
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

alpha0 = 2.0;
coefficientBound = 50;

z0 = [0.9; -0.6];
Tengage = 1.0;
Tfinal = 12.0;

saveFigures = true;

%% Model parameters

par.Je = 1.0;
par.Jd = 1.5;

% Mode 1: disengaged/open-loop dynamics.
par.a1_open = 0.20;
par.b1_open = 0.15;
par.a2_open = 0.10;
par.b2_open = 0.05;

% Mode 2: slipping dynamics.
par.a1 = 0.20;
par.b1 = 0.10;
par.a2 = 0.20;
par.b2 = 0.05;
par.c1 = 0.20;
par.c3 = 0.03;
par.gamma = par.Je/par.Jd;

% Mode 3: locked scalar dynamics.
par.al  = 0.30;
par.bl  = 0.05;
par.eta = 1/(par.Je + par.Jd);

% Locking threshold.
par.delta = 0.08;

%% Solver options

opts = sdpsettings( ...
    'solver','mosek', ...
    'verbose',1, ...
    'sos.model',2);

%% Book plotting style

S.blue  = [0, 92, 175]/255;
S.black = [0, 0, 0];

%% Symbolic variables

sdpvar z1 z2 y
z = [z1; z2];

%% Fixed Lyapunov functions

Vslip = 0.5*(par.Je*z1^2 + par.Jd*z2^2);
Vlock = 0.5*(par.Je + par.Jd)*y^2;

%% Linear feedback laws

sdpvar ae1 ae2 ac1 ac2 alock

ke = ae1*z1 + ae2*z2;
kc = ac1*z1 + ac2*z2;
klock = alock*y;

%% Slipping-mode dynamics

slip = z1 - z2;

tauPassive = par.c1*slip + par.c3*slip^3;
tauTotal = tauPassive + kc;

fSlip = [
    par.a1*z1 - par.b1*z1^3 + ke - tauTotal;
   -par.a2*z2 - par.b2*z2^3 + par.gamma*tauTotal
];

LVslip = jacobian(Vslip,z)*fSlip;

%% Locked-mode dynamics

fLock = -par.al*y - par.bl*y^3 + par.eta*klock;
LVlock = jacobian(Vlock,y)*fLock;

%% Reset compatibility: Vlock(R(z)) <= Vslip(z)

zbar = (par.Je*z1 + par.Jd*z2)/(par.Je + par.Jd);

VlockReset = 0.5*(par.Je + par.Jd)*zbar^2;
jumpPolynomial = clean(Vslip - VlockReset,1e-12);

%% SOS synthesis problem

constraints = [];

constraints = [constraints, sos(-LVslip - alpha0*Vslip)];
constraints = [constraints, sos(-LVlock - alpha0*Vlock)];
constraints = [constraints, sos(jumpPolynomial)];

constraints = [constraints, -coefficientBound <= ae1 <= coefficientBound];
constraints = [constraints, -coefficientBound <= ae2 <= coefficientBound];
constraints = [constraints, -coefficientBound <= ac1 <= coefficientBound];
constraints = [constraints, -coefficientBound <= ac2 <= coefficientBound];
constraints = [constraints, -coefficientBound <= alock <= coefficientBound];

objective = ae1^2 + ae2^2 + ac1^2 + ac2^2 + alock^2;

fprintf('\n============================================================\n');
fprintf('Hybrid clutch SOS controller synthesis\n');
fprintf('============================================================\n\n');

fprintf('Prescribed decay rate alpha0 = %.6g\n',alpha0);
fprintf('Solving minimum-gain SOS synthesis problem...\n');

sol = solvesos(constraints,objective,opts);

if sol.problem ~= 0
    error('SOS synthesis failed: %s',sol.info);
end

ae1Value = value(ae1);
ae2Value = value(ae2);
ac1Value = value(ac1);
ac2Value = value(ac2);
alockValue = value(alock);

fprintf('\nSynthesized controllers:\n');
fprintf('  ke(z1,z2) = %.12g*z1 %+ .12g*z2\n',ae1Value,ae2Value);
fprintf('  kc(z1,z2) = %.12g*z1 %+ .12g*z2\n',ac1Value,ac2Value);
fprintf('  k3(y)     = %.12g*y\n',alockValue);

fprintf('\nJump polynomial Vslip - Vlock(R(z)):\n');
fprintf('  coefficient of (z1-z2)^2 = %.12g\n', ...
    par.Je*par.Jd/(2*(par.Je+par.Jd)));

keFun = @(z1,z2) ae1Value*z1 + ae2Value*z2;
kcFun = @(z1,z2) ac1Value*z1 + ac2Value*z2;
klockFun = @(y) alockValue*y;

%% Simulate uncontrolled and controlled hybrid trajectories

[tUncontrolled,zUncontrolled,modeUncontrolled] = simulate_hybrid_clutch( ...
    z0,Tengage,Tfinal,par,[],[],[],false);

[tControlled,zControlled,modeControlled] = simulate_hybrid_clutch( ...
    z0,Tengage,Tfinal,par,keFun,kcFun,klockFun,true);

%% Figure 1: slip synchronization

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 12 8]);

hold on
box on

plot(tUncontrolled,zUncontrolled(:,1)-zUncontrolled(:,2), ...
    '--', ...
    'Color',S.black, ...
    'LineWidth',2);

plot(tControlled,zControlled(:,1)-zControlled(:,2), ...
    '-', ...
    'Color',S.blue, ...
    'LineWidth',4);

yline(par.delta,':k','LineWidth',2);
yline(-par.delta,':k','LineWidth',2);

xlabel('$t$','Interpreter','latex','FontSize',36);
ylabel('$z_1-z_2$','Interpreter','latex','FontSize',36);

legend({'uncontrolled','controlled','$\pm\delta$'}, ...
    'Interpreter','latex', ...
    'Location','best');

set(gca,'FontSize',24,'TickLabelInterpreter','latex','LineWidth',1.0);

xlim([0 Tfinal]);
ylim([-0.15 1.5]);

grid on

if saveFigures
    export_pdf(fig1,'hybrid_control_1.pdf');
end

%% Figure 2: hybrid mode evolution

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 12 8]);

hold on
box on

stairs(tUncontrolled,modeUncontrolled, ...
    '--', ...
    'Color',S.black, ...
    'LineWidth',2);

stairs(tControlled,modeControlled, ...
    '-', ...
    'Color',S.blue, ...
    'LineWidth',4);

yticks([1 2 3]);
yticklabels({'open','slip','locked'});

xlabel('$t$','Interpreter','latex','FontSize',36);
ylabel('mode','Interpreter','latex','FontSize',36);

legend('uncontrolled','controlled','Location','best');

set(gca,'FontSize',24,'TickLabelInterpreter','latex','LineWidth',1.0);

xlim([0 Tfinal]);

grid on

if saveFigures
    export_pdf(fig2,'hybrid_control_2.pdf');
end

fprintf('\nFinished hybrid clutch control example.\n');

%% Local functions

function [tout,zout,modeout] = simulate_hybrid_clutch( ...
    z0,Tengage,Tfinal,par,keFun,kcFun,klockFun,useControl)

    tout = [];
    zout = [];
    modeout = [];

    zcur = z0(:);

    %% Mode 1: open

    [t,z] = ode45( ...
        @(t,z) hybrid_rhs(1,z,par,keFun,kcFun,klockFun,useControl), ...
        [0 Tengage], ...
        zcur);

    append_segment(1,t,z);

    zcur = z(end,:).';

    %% Mode 2: slipping

    opts = odeset('Events',@(t,z) lock_event(t,z,par));

    [t,z] = ode45( ...
        @(t,z) hybrid_rhs(2,z,par,keFun,kcFun,klockFun,useControl), ...
        [Tengage Tfinal], ...
        zcur, ...
        opts);

    append_segment(2,t,z);

    zcur = z(end,:).';
    t0 = t(end);

    %% Mode 3: locked

    if t0 < Tfinal

        zbar = (par.Je*zcur(1) + par.Jd*zcur(2))/(par.Je+par.Jd);
        zcur = [zbar; zbar];

        [t,z] = ode45( ...
            @(t,z) hybrid_rhs(3,z,par,keFun,kcFun,klockFun,useControl), ...
            [t0 Tfinal], ...
            zcur);

        append_segment(3,t,z);
    end

    function append_segment(mode,t,z)

        if isempty(tout)

            tout = t;
            zout = z;
            modeout = mode*ones(size(t));

        else

            tout = [tout; t(2:end)]; 
            zout = [zout; z(2:end,:)]; 
            modeout = [modeout; mode*ones(length(t)-1,1)]; 
        end
    end
end

function dz = hybrid_rhs(mode,z,par,keFun,kcFun,klockFun,useControl)

    z1 = z(1);
    z2 = z(2);

    switch mode

        case 1

            dz = [
                par.a1_open*z1 - par.b1_open*z1^3;
               -par.a2_open*z2 - par.b2_open*z2^3
            ];

        case 2

            if useControl
                ue = keFun(z1,z2);
                uc = kcFun(z1,z2);
            else
                ue = 0;
                uc = 0;
            end

            slip = z1 - z2;

            tauPassive = par.c1*slip + par.c3*slip^3;
            tauTotal = tauPassive + uc;

            dz = [
                par.a1*z1 - par.b1*z1^3 + ue - tauTotal;
               -par.a2*z2 - par.b2*z2^3 + par.gamma*tauTotal
            ];

        case 3

            y = z1;

            if useControl
                u = klockFun(y);
            else
                u = 0;
            end

            dy = -par.al*y - par.bl*y^3 + par.eta*u;

            dz = [dy; dy];
    end
end

function [value,isterminal,direction] = lock_event(~,z,par)

    value = (z(1)-z(2))^2 - par.delta^2;
    isterminal = 1;
    direction = -1;
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