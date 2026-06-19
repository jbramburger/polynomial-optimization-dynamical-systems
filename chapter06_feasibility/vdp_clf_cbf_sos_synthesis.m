%% vdp_clf_cbf_sos_synthesis.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: CLF-CBF synthesis for the controlled Van der Pol oscillator
%
% This script synthesizes polynomial feedback controllers for the controlled
% Van der Pol oscillator
%
%     x1' = x2,
%     x2' = mu (1 - x1^2) x2 - x1 + u.
%
% A quadratic control Lyapunov function V is prescribed. For each certified
% safe set, the script searches for a polynomial feedback law u = k(x) such
% that:
%
%   1. V decreases on a compact certification region, and
%   2. the safe set is forward invariant by control-barrier conditions.
%
% Two safe sets are considered:
%
%   1. a velocity-constrained ellipse,
%   2. a three-lobed semialgebraic region.
%
% Pressing Run should reproduce the numerical output and save the figures
%
%     vdp_clf_cbf_ellipse.pdf
%     vdp_clf_cbf_lobe.pdf
%
% in the current directory.
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

mu = 1;

eps_clf = 1e-4;
beta_cbf = 1.0;

RK = 2.6;

degK = 3;

degSigCLF   = 4;
degSigK_CLF = 2;

degSigCBF   = 4;
degSigK_CBF = 2;
degSigH_CBF = 2;

solver_name = 'mosek';
verbose_solver = 1;

saveFigures = true;

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.black  = [0, 0, 0];
S.gray   = [0.40, 0.40, 0.40];

%% Symbolic variables and system

sdpvar x1 x2
x = [x1; x2];

f = [
    x2;
    mu*(1 - x1^2)*x2 - x1
];

g = [0; 1];

V = 1.5*x1^2 + x1*x2 + 0.5*x2^2;
p = x1^2 + x2^2;

qK = RK^2 - x1^2 - x2^2;

fprintf('\n============================================================\n');
fprintf('CLF-CBF synthesis for the controlled Van der Pol oscillator\n');
fprintf('============================================================\n\n');

fprintf('mu       = %.6g\n', mu);
fprintf('RK       = %.6g\n', RK);
fprintf('degK     = %d\n', degK);
fprintf('beta_cbf = %.6g\n\n', beta_cbf);

%% Barrier 1: velocity-constrained ellipse

a = 2.0;
b = 0.8;

hEllipse = 1 - x1^2/a^2 - x2^2/b^2;

[kEllipse, successEllipse] = synthesize_controller( ...
    x, f, g, V, p, qK, hEllipse, eps_clf, beta_cbf, ...
    degK, degSigCLF, degSigK_CLF, degSigCBF, degSigK_CBF, degSigH_CBF, ...
    solver_name, verbose_solver, 'ellipse');

if ~successEllipse
    error('Ellipse barrier controller synthesis failed.');
end

%% Barrier 2: three-lobed semialgebraic safe region

eta = 0.1;
cLobe = 3.0;

hLobe = cLobe - V + eta*(x1^3 - 3*x1*x2^2);

[kLobe, successLobe] = synthesize_controller( ...
    x, f, g, V, p, qK, hLobe, eps_clf, beta_cbf, ...
    degK, degSigCLF, degSigK_CLF, degSigCBF, degSigK_CBF, degSigH_CBF, ...
    solver_name, verbose_solver, 'three-lobed');

if ~successLobe
    error('Three-lobed barrier controller synthesis failed.');
end

%% Display controllers

fprintf('\nSynthesized controller for ellipse barrier:\n');
disp(sdisplay(clean(kEllipse,1e-7)))

fprintf('\nSynthesized controller for three-lobed barrier:\n');
disp(sdisplay(clean(kLobe,1e-7)))

%% Plot settings

xmin = -2.6; xmax = 2.6;
ymin = -2.6; ymax = 2.6;

nx = 300;
ny = 300;

[xg, yg] = meshgrid(linspace(xmin,xmax,nx), linspace(ymin,ymax,ny));

VV = 1.5*xg.^2 + xg.*yg + 0.5*yg.^2;

HEllipse = 1 - xg.^2/a^2 - yg.^2/b^2;
HLobe = cLobe - VV + eta*(xg.^3 - 3*xg.*yg.^2);
QKgrid = RK^2 - xg.^2 - yg.^2;

%% Initial conditions

ICellipse = [
    -1.55  0.18
    -1.20 -0.35
    -0.70  0.55
     0.70 -0.55
     1.20  0.35
     1.55 -0.18
];

h0e = 1 - ICellipse(:,1).^2/a^2 - ICellipse(:,2).^2/b^2;
q0e = RK^2 - ICellipse(:,1).^2 - ICellipse(:,2).^2;

if any(h0e <= 0) || any(q0e <= 0)
    error('At least one ellipse initial condition is outside S_cert.');
end

IClobe = [
     0.00 -2.00
    -0.80  2.20
    -0.95 -0.95
    -1.40  0.50
    -1.70  1.60
     0.82 -2.20
     0.10  1.90
     1.30 -0.65
];

V0l = 1.5*IClobe(:,1).^2 ...
    + IClobe(:,1).*IClobe(:,2) ...
    + 0.5*IClobe(:,2).^2;

h0l = cLobe - V0l ...
    + eta*(IClobe(:,1).^3 - 3*IClobe(:,1).*IClobe(:,2).^2);

q0l = RK^2 - IClobe(:,1).^2 - IClobe(:,2).^2;

if any(h0l <= 0) || any(q0l <= 0)
    error('At least one lobe initial condition is outside S_cert.');
end

%% Simulate trajectories

odeopts = odeset('RelTol',1e-9,'AbsTol',1e-11);
Tfwd = [0 15];

trajEllipse = simulate_trajectories(ICellipse,kEllipse,x,mu,Tfwd,odeopts);
trajLobe    = simulate_trajectories(IClobe,kLobe,x,mu,Tfwd,odeopts);

%% Figure A: ellipse safe set

figA = figure;
set(figA,'Color','w','Units','centimeters','Position',[2 2 12 11]);

hold on
box on

contourf(xg,yg,VV,18,'LineStyle','none');
colormap(gca,parula);

contour(xg,yg,VV,14, ...
    'Color',[0.85 0.85 0.85], ...
    'LineWidth',0.55);

contourf(xg,yg,double(HEllipse >= 0 & QKgrid >= 0),[0.5 0.5], ...
    'LineStyle','none', ...
    'FaceAlpha',0.10);

contour(xg,yg,HEllipse,[0 0], ...
    'Color',S.black, ...
    'LineWidth',5);

th = linspace(0,2*pi,500);

plot(RK*cos(th),RK*sin(th),'--', ...
    'Color',S.gray, ...
    'LineWidth',5);

for j = 1:length(trajEllipse)

    ztraj = trajEllipse{j};

    plot(ztraj(:,1),ztraj(:,2),'-', ...
        'Color',S.red, ...
        'LineWidth',5);

    plot(ztraj(1,1),ztraj(1,2),'o', ...
        'MarkerSize',5.2, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',S.red, ...
        'LineWidth',5);
end

plot(0,0,'o', ...
    'MarkerSize',20, ...
    'MarkerFaceColor',S.green, ...
    'MarkerEdgeColor',S.black, ...
    'LineWidth',3);

axis([xmin xmax ymin ymax]);

set(gca,'FontSize',42)
xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

if saveFigures
    save_pdf_figure(figA,'vdp_clf_cbf_ellipse.pdf');
end

%% Figure B: three-lobed safe set

figB = figure;
set(figB,'Color','w','Units','centimeters','Position',[2 2 12 11]);

hold on
box on

contourf(xg,yg,VV,18,'LineStyle','none');
colormap(gca,parula);

contour(xg,yg,VV,14, ...
    'Color',[0.85 0.85 0.85], ...
    'LineWidth',0.55);

contourf(xg,yg,double(HLobe >= 0 & QKgrid >= 0),[0.5 0.5], ...
    'LineStyle','none', ...
    'FaceAlpha',0.10);

contour(xg,yg,HLobe,[0 0], ...
    'Color',S.black, ...
    'LineWidth',5);

plot(RK*cos(th),RK*sin(th),'--', ...
    'Color',S.gray, ...
    'LineWidth',5);

for j = 1:length(trajLobe)

    ztraj = trajLobe{j};

    plot(ztraj(:,1),ztraj(:,2),'-', ...
        'Color',S.red, ...
        'LineWidth',5);

    plot(ztraj(1,1),ztraj(1,2),'o', ...
        'MarkerSize',5.2, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',S.red, ...
        'LineWidth',5);
end

plot(0,0,'o', ...
    'MarkerSize',20, ...
    'MarkerFaceColor',S.green, ...
    'MarkerEdgeColor',S.black, ...
    'LineWidth',3);

axis([xmin xmax ymin ymax]);

set(gca,'FontSize',42)
xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

if saveFigures
    save_pdf_figure(figB,'vdp_clf_cbf_lobe.pdf');
end

fprintf('\nExported figures:\n');
fprintf('  vdp_clf_cbf_ellipse.pdf\n');
fprintf('  vdp_clf_cbf_lobe.pdf\n');
fprintf('\nFinished Van der Pol CLF-CBF synthesis example.\n');

%% Local functions

function [kfinal,success] = synthesize_controller( ...
    x, f, g, V, p, qK, h, eps_clf, beta_cbf, ...
    degK, degSigCLF, degSigK_CLF, degSigCBF, degSigK_CBF, degSigH_CBF, ...
    solver_name, verbose_solver, tag)

    fprintf('\nSolving SOS feasibility problem for %s barrier...\n',tag);

    x1 = x(1);
    x2 = x(2);

    [k,cK] = polynomial(x,degK);

    k0 = replace(k,x,[0;0]);
    fcl = f + g*k;

    Vdot = jacobian(V,x1)*fcl(1) + jacobian(V,x2)*fcl(2);
    hdot = jacobian(h,x1)*fcl(1) + jacobian(h,x2)*fcl(2);
    qdot = jacobian(qK,x1)*fcl(1) + jacobian(qK,x2)*fcl(2);

    [sigCLF,cSigCLF] = polynomial(x,degSigCLF);
    [sigK_CLF,cSigK_CLF] = polynomial(x,degSigK_CLF);

    rCLF = -Vdot - eps_clf*p - sigCLF - sigK_CLF*qK;

    [sigH0,cSigH0] = polynomial(x,degSigCBF);
    [sigHK,cSigHK] = polynomial(x,degSigK_CBF);
    [sigHH,cSigHH] = polynomial(x,degSigH_CBF);

    rH = hdot + beta_cbf*h ...
        - sigH0 - sigHK*qK - sigHH*h;

    [sigQ0,cSigQ0] = polynomial(x,degSigCBF);
    [sigQK,cSigQK] = polynomial(x,degSigK_CBF);
    [sigQH,cSigQH] = polynomial(x,degSigH_CBF);

    rQ = qdot + beta_cbf*qK ...
        - sigQ0 - sigQK*qK - sigQH*h;

    constraints = [
        k0 == 0, ...
        sos(sigCLF), ...
        sos(sigK_CLF), ...
        sos(sigH0), ...
        sos(sigHK), ...
        sos(sigHH), ...
        sos(sigQ0), ...
        sos(sigQK), ...
        sos(sigQH), ...
        coefficients(rCLF,x) == 0, ...
        coefficients(rH,x) == 0, ...
        coefficients(rQ,x) == 0
        ];

    opts = sdpsettings('solver',solver_name,'verbose',verbose_solver);

    sol = solvesos(constraints, [], opts, ...
        [cK; ...
         cSigCLF; cSigK_CLF; ...
         cSigH0; cSigHK; cSigHH; ...
         cSigQ0; cSigQK; cSigQH]);

    success = (sol.problem == 0);

    if ~success
        fprintf('SOS feasibility failed for %s barrier: %s\n',tag,sol.info);
        kfinal = NaN;
        return;
    end

    cKval = value(cK);

    if any(isnan(cKval))
        fprintf('NaN controller coefficients returned for %s barrier.\n',tag);
        success = false;
        kfinal = NaN;
        return;
    end

    kfinal = clean(replace(k,cK,cKval),1e-8);

    fprintf('SOS feasibility succeeded for %s barrier.\n',tag);
end

function traj = simulate_trajectories(IC,kpoly,x,mu,Tfwd,opts)

    traj = cell(size(IC,1),1);

    for j = 1:size(IC,1)

        z0 = IC(j,:)';
        ode = @(t,z) controlled_vdp_rhs(t,z,kpoly,x,mu);

        [~,zz] = ode45(ode,Tfwd,z0,opts);
        traj{j} = zz;
    end
end

function dz = controlled_vdp_rhs(~,z,kpoly,x,mu)

    x1 = z(1);
    x2 = z(2);

    u = eval_yalmip_poly_point(kpoly,x,x1,x2);

    dz = [
        x2;
        mu*(1 - x1^2)*x2 - x1 + u
    ];
end

function val = eval_yalmip_poly_point(poly,x,a,b)

    [coef,mon] = coefficients(poly,x);
    coef = value(coef);

    val = 0;

    for kk = 1:length(coef)

        if abs(coef(kk)) > 1e-12

            e1 = degree(mon(kk),x(1));
            e2 = degree(mon(kk),x(2));

            val = val + coef(kk)*(a^e1)*(b^e2);
        end
    end
end

function save_pdf_figure(fig,filename)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',filename);
    fprintf('Saved figure: %s\n', filename);
end