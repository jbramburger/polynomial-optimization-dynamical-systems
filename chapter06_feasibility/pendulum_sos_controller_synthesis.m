%% pendulum_sos_controller_synthesis.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: SOS controller synthesis for the upright torque-actuated pendulum
%
% This script synthesizes a polynomial feedback controller for the upright
% torque-actuated pendulum
%
%     theta' = omega,
%     omega' = -sin(theta) + u.
%
% The dynamics are lifted using
%
%     x1 = sin(theta),   x2 = cos(theta),   x3 = omega,
%
% and shifted so that the upright equilibrium theta = pi, omega = 0 becomes
% the origin in the variables
%
%     z1 = x1,   z2 = x2 + 1,   z3 = x3.
%
% The physical manifold is
%
%     h(z) = 2 z2 - z1^2 - z2^2 = 0.
%
% A quadratic Lyapunov function V is prescribed, and the script searches for
% a polynomial feedback law u = k(z) such that V decreases on the physical
% manifold inside a certified local sublevel set.
%
% Pressing Run should reproduce the numerical output and save the figures
%
%     pendulum_panelA.pdf
%     pendulum_panelB.pdf
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

eps_clf = 1e-5;

degK = 3;
degSigmaDec = 4;
degSigmaQDec = 2;
degLambda = 2;

cLevelCert = 3;

solver_name = 'mosek';
verbose_solver = 1;

saveFigures = true;

opts = sdpsettings('solver',solver_name,'verbose',verbose_solver);

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.black  = [0, 0, 0];

%% YALMIP variables

sdpvar z1 z2 z3
z = [z1; z2; z3];

h = 2*z2 - z1^2 - z2^2;

p = z1^2 + z2^2 + z3^2;

%% Prescribed Lyapunov function

Vfinal = z1^2 + z2^2 + z3^2 - z1*z3;

q = cLevelCert - Vfinal;

%% Polynomial controller decision variable

[k,cK] = polynomial(z,degK);

k0 = replace(k,z,[0;0;0]);

%% Closed-loop lifted dynamics

fcl = [
    (z2 - 1)*z3;
    -z1*z3;
    -z1 + k
];

%% SOS decay certificate

[sig1,cSig1] = polynomial(z,degSigmaDec);
[sigq1,cSigq1] = polynomial(z,degSigmaQDec);
[lam1,cLam1] = polynomial(z,degLambda);

Vdot = jacobian(Vfinal,z)*fcl;

rDec = -Vdot - eps_clf*p - sig1 - sigq1*q - lam1*h;

constraints = [
    k0 == 0, ...
    sos(sig1), ...
    sos(sigq1), ...
    coefficients(rDec,z) == 0
    ];

objective = 1e-6*(cK'*cK + cLam1'*cLam1);

fprintf('\n============================================================\n');
fprintf('SOS controller synthesis for the upright pendulum\n');
fprintf('============================================================\n\n');

sol = solvesos(constraints, objective, opts, ...
    [cK; cSig1; cSigq1; cLam1]);

if sol.problem ~= 0
    error('SOS controller synthesis failed: %s', sol.info);
end

fprintf('SOS controller synthesis successful.\n');

cKval = value(cK);

if any(isnan(cKval))
    error('Some coefficients of k are NaN. Increase regularization or constrain the ansatz.');
end

kfinal = replace(k,cK,cKval);
kfinal = clean(kfinal,1e-8);

disp('Prescribed Lyapunov function V(z):');
disp(sdisplay(Vfinal))

disp('Synthesized feedback law k(z):');
disp(sdisplay(kfinal))

%% Extract approximate linear gains

k_z1 = value(replace(jacobian(kfinal,z1),z,[0;0;0]));
k_z2 = value(replace(jacobian(kfinal,z2),z,[0;0;0]));
k_z3 = value(replace(jacobian(kfinal,z3),z,[0;0;0]));

fprintf('\nLinearization of synthesized controller near upright:\n');
fprintf('u ≈ %.8f z1 + %.8f z2 + %.8f z3\n', k_z1, k_z2, k_z3);

%% Grid in original pendulum coordinates

phiGrid = linspace(-2.2,2.2,300);
omGrid  = linspace(-2.2,2.2,220);

[PHI,OMG] = meshgrid(phiGrid,omGrid);

TH = pi + PHI;

Z1 = sin(TH);
Z2 = cos(TH) + 1;
Z3 = OMG;

%% Evaluate V and Vdot on the grid

VV = eval_yalmip_poly_on_grid(Vfinal,z,Z1,Z2,Z3);

dVdz = jacobian(Vfinal,z);

VX = eval_yalmip_poly_on_grid(dVdz(1),z,Z1,Z2,Z3);
VY = eval_yalmip_poly_on_grid(dVdz(2),z,Z1,Z2,Z3);
VZ = eval_yalmip_poly_on_grid(dVdz(3),z,Z1,Z2,Z3);

U = eval_yalmip_poly_on_grid(kfinal,z,Z1,Z2,Z3);

F1 = (Z2 - 1).*Z3;
F2 = -Z1.*Z3;
F3 = -Z1 + U;

VDOT = VX.*F1 + VY.*F2 + VZ.*F3; 

%% Simulate trajectories

odefun = @(t,Y) pendulum_closed_loop_synth(t,Y,kfinal,z);

phi0 = [-1.2 -0.9 -0.55 0.55 0.9 1.2];
omega0 = [0.6 -1.0 1.3 -1.3 1.0 -0.6];

theta0 = pi + phi0;

z1_0 = sin(theta0);
z2_0 = cos(theta0) + 1;
z3_0 = omega0;

V0 = z1_0.^2 + z2_0.^2 + z3_0.^2 - z1_0.*z3_0;

if any(V0 >= cLevelCert)
    error('At least one initial condition lies outside V <= cLevelCert.');
end

fprintf('\nMaximum initial V-value: %.6f\n', max(V0));

T = [0 15];

traj = cell(length(theta0),1);

for j = 1:length(theta0)
    [~,YY] = ode45(odefun,T,[theta0(j);omega0(j)]);
    traj{j} = YY;
end

%% Figure A: contours of V and closed-loop trajectories

figA = figure;
set(figA,'Color','w','Units','centimeters','Position',[2 2 12 11]);

hold on
box on

contourf(PHI,OMG,VV,18,'LineStyle','none');
colormap(gca, parula);

contour(PHI,OMG,VV,14, ...
    'Color',[0.85 0.85 0.85], ...
    'LineWidth',0.55);

contour(PHI,OMG,VV,[cLevelCert cLevelCert], ...
    'Color',S.black, ...
    'LineWidth',3.0);

contourf(PHI,OMG,VV,[min(VV(:)) cLevelCert], ...
    'LineStyle','none', ...
    'FaceAlpha',0.15);

for j = 1:length(traj)

    YY = traj{j};

    plot(wrapToPiLocal(YY(:,1)-pi),YY(:,2), ...
        'Color',S.red, ...
        'LineWidth',7);

    plot(wrapToPiLocal(YY(1,1)-pi),YY(1,2), ...
        'o', ...
        'MarkerSize',5.8, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',S.red, ...
        'LineWidth',5);
end

plot(0,0,'o', ...
    'MarkerSize',20, ...
    'MarkerFaceColor',S.green, ...
    'MarkerEdgeColor',S.black, ...
    'LineWidth',1.0);

axis([-2.2 2.2 -2.2 2.2]);

set(gca,'FontSize',42)
xlabel('$\theta-\pi$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$\omega$','Interpreter','latex','FontSize',42,'FontWeight','bold');

colorbar;

if saveFigures
    save_pdf_figure(figA,'pendulum_panelA.pdf');
end

%% Figure B: 3D Lyapunov landscape

figB = figure;
set(figB,'Color','w','Units','centimeters','Position',[2 2 12 11]);

hold on
box on
grid on

Vcap = prctile(VV(:),80);
VVplot = VV;
VVplot(VVplot > Vcap) = Vcap;

surf(PHI,OMG,VVplot, ...
    'EdgeColor','none', ...
    'FaceAlpha',1);

colormap(parula);

C = contourc(phiGrid,omGrid,VV,[cLevelCert cLevelCert]);

idx = 1;
while idx < size(C,2)

    npts = C(2,idx);
    pts = C(:,idx+1:idx+npts);

    phiC = pts(1,:);
    omC  = pts(2,:);

    plot3(phiC,omC,cLevelCert*ones(size(phiC)), ...
        'Color',S.black, ...
        'LineWidth',4.0);

    idx = idx+npts+1;
end

plot3(0,0,0,'o', ...
    'MarkerSize',10, ...
    'MarkerFaceColor',S.green, ...
    'MarkerEdgeColor',S.black, ...
    'LineWidth',1.0);

axis([-2.2 2.2 -2.2 2.2 0 Vcap]);
view(42,28);

set(gca,'FontSize',42)
xlabel('$\theta-\pi$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$\omega$','Interpreter','latex','FontSize',42,'FontWeight','bold');
zlabel('$V$','Interpreter','latex','FontSize',42,'FontWeight','bold');

if saveFigures
    save_pdf_figure(figB,'pendulum_panelB.pdf');
end

fprintf('\nFinished pendulum SOS controller synthesis example.\n');

%% Local functions

function dY = pendulum_closed_loop_synth(~,Y,kpoly,z)

    theta = Y(1);
    omega = Y(2);

    zz1 = sin(theta);
    zz2 = cos(theta) + 1;
    zz3 = omega;

    u = eval_yalmip_poly_point(kpoly,z,zz1,zz2,zz3);

    dtheta = omega;
    domega = -sin(theta) + u;

    dY = [dtheta; domega];
end

function y = wrapToPiLocal(x)

    y = mod(x + pi,2*pi) - pi;
end

function F = eval_yalmip_poly_on_grid(poly,z,Z1,Z2,Z3)

    [coef,mon] = coefficients(poly,z);
    coef = value(coef);

    F = zeros(size(Z1));

    for kk = 1:length(coef)

        if abs(coef(kk)) > 1e-12

            e1 = degree(mon(kk),z(1));
            e2 = degree(mon(kk),z(2));
            e3 = degree(mon(kk),z(3));

            F = F + coef(kk).*(Z1.^e1).*(Z2.^e2).*(Z3.^e3);
        end
    end
end

function val = eval_yalmip_poly_point(poly,z,a,b,c0)

    [coef,mon] = coefficients(poly,z);
    coef = value(coef);

    val = 0;

    for kk = 1:length(coef)

        if abs(coef(kk)) > 1e-12

            e1 = degree(mon(kk),z(1));
            e2 = degree(mon(kk),z(2));
            e3 = degree(mon(kk),z(3));

            val = val + coef(kk)*(a^e1)*(b^e2)*(c0^e3);
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