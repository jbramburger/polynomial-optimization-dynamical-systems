%% sos_lyapunov_homogeneous_cubic_full.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: SOS Lyapunov functions for a homogeneous cubic system
%
% This script studies the planar homogeneous cubic system
%
%     x' = -sin(theta) x^3 + cos(theta) y^3,
%     y' = -cos(theta) x^3 - sin(theta) y^3.
%
% The script:
%
%   1. Displays the exact quartic Lyapunov function V = x^4 + y^4.
%   2. Checks infeasibility of a quadratic SOS Lyapunov search.
%   3. Computes a quartic SOS Lyapunov function.
%   4. Checks a homogeneous sextic SOS Lyapunov search.
%   5. Computes a nonhomogeneous sextic SOS Lyapunov function.
%   6. Produces the figures used in the book.
%
% Pressing Run should reproduce the numerical output and save the figures
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

theta       = 0.01;
eps_sos     = 1e-3;
solver_name = 'mosek';
saveFigures = true;

xmin = -1.8; xmax = 1.8;
ymin = -1.8; ymax = 1.8;
Nx = 401; Ny = 401;

opts = sdpsettings('solver',solver_name,'verbose',0);

%% Book plotting style

blue   = [0, 92, 175]/255;
red    = [200, 50, 50]/255;
green  = [0, 140, 90]/255;
purple = [170, 90, 160]/255;
black  = [0, 0, 0];

lw_main = 2.4;

%% Variables and vector field

sdpvar x y

s = sin(theta);
c = cos(theta);

f1 = -s*x^3 + c*y^3;
f2 = -c*x^3 - s*y^3;
f  = [f1; f2];

r2 = x^2 + y^2;

%% Exact quartic Lyapunov function

V_exact = x^4 + y^4;
Vdot_exact = clean(jacobian(V_exact,[x y])*f,1e-12);

fprintf('\n============================================================\n');
fprintf('SOS Lyapunov functions for a homogeneous cubic system\n');
fprintf('============================================================\n\n');

fprintf('Parameter values:\n');
fprintf('  theta   = %.12g\n', theta);
fprintf('  eps_sos = %.12g\n\n', eps_sos);

fprintf('Exact quartic Lyapunov function:\n');
fprintf('  V_exact(x,y) = x^4 + y^4\n');
fprintf('  Vdot_exact   = -4 sin(theta) (x^6 + y^6)\n');
fprintf('               = %.12g*x^6 + %.12g*y^6\n\n', ...
    -4*sin(theta), -4*sin(theta));

%% Quadratic SOS Lyapunov search

fprintf('------------------------------------------------------------\n');
fprintf('Checking quadratic SOS Lyapunov search\n');
fprintf('------------------------------------------------------------\n');

sdpvar a2 b2 c2

W2 = a2*x^2 + b2*x*y + c2*y^2;
W2dot = jacobian(W2,[x y])*f;

F2 = [
    sos(W2 - eps_sos*r2), ...
    sos(-W2dot - eps_sos*r2^2), ...
    replace(W2,[x y],[1 0]) == 1
    ];

sol2 = solvesos(F2,[],opts,[a2;b2;c2]);

fprintf('Quadratic status code: %d\n', sol2.problem);
fprintf('Quadratic solver info: %s\n', sol2.info);

quadratic_feasible = (sol2.problem == 0);

if quadratic_feasible
    fprintf('Quadratic SOS search feasible.\n');
    fprintf('Recovered quadratic:\n');
    fprintf('  W2(x,y) = %.12g*x^2 + %.12g*x*y + %.12g*y^2\n\n', ...
        value(a2), value(b2), value(c2));
else
    fprintf('Quadratic SOS search infeasible or numerically rejected.\n\n');
end

%% Quartic SOS Lyapunov search

fprintf('------------------------------------------------------------\n');
fprintf('Searching for quartic SOS Lyapunov function\n');
fprintf('------------------------------------------------------------\n');

sdpvar a4 b4

W4 = a4*(x^4 + y^4) + b4*x^2*y^2;
W4dot = jacobian(W4,[x y])*f;

F4 = [
    sos(W4 - eps_sos*r2^2), ...
    sos(-W4dot - eps_sos*r2^3), ...
    replace(W4,[x y],[1 0]) == 1
    ];

sol4 = solvesos(F4,[],opts,[a4;b4]);

fprintf('Quartic status code: %d\n', sol4.problem);
fprintf('Quartic solver info: %s\n', sol4.info);

quartic_feasible = (sol4.problem == 0);

if ~quartic_feasible
    error('Quartic SOS search failed. Solver message: %s', sol4.info);
end

a4sol = value(a4);
b4sol = value(b4);

fprintf('Quartic SOS search feasible.\n');
fprintf('Recovered quartic Lyapunov function:\n');
fprintf('  W4(x,y) = %.12g*(x^4+y^4) + %.12g*x^2*y^2\n\n', ...
    a4sol, b4sol);

%% Homogeneous sextic SOS Lyapunov search

fprintf('------------------------------------------------------------\n');
fprintf('Checking homogeneous sextic SOS Lyapunov search\n');
fprintf('------------------------------------------------------------\n');

sdpvar ah6 bh6

H6 = ah6*(x^6 + y^6) + bh6*(x^4*y^2 + x^2*y^4);
H6dot = jacobian(H6,[x y])*f;

FH6 = [
    sos(H6 - eps_sos*r2^3), ...
    sos(-H6dot - eps_sos*r2^4), ...
    replace(H6,[x y],[1 0]) == 1
    ];

solH6 = solvesos(FH6,[],opts,[ah6;bh6]);

fprintf('Homogeneous sextic status code: %d\n', solH6.problem);
fprintf('Homogeneous sextic solver info: %s\n', solH6.info);

homogeneous_sextic_feasible = (solH6.problem == 0);

if homogeneous_sextic_feasible
    ah6sol = value(ah6);
    bh6sol = value(bh6);

    fprintf('Homogeneous sextic SOS search feasible.\n');
    fprintf('Recovered homogeneous sextic:\n');
    fprintf('  H6(x,y) = %.12g*(x^6+y^6) + %.12g*(x^4*y^2+x^2*y^4)\n\n', ...
        ah6sol, bh6sol);
else
    fprintf('Homogeneous sextic SOS search infeasible or numerically rejected.\n\n');
end

%% Nonhomogeneous sextic SOS Lyapunov search

fprintf('------------------------------------------------------------\n');
fprintf('Searching for nonhomogeneous sextic SOS Lyapunov function\n');
fprintf('------------------------------------------------------------\n');

sdpvar a6 b6 c6 d6 e6

W6 = a6*(x^6 + y^6) ...
   + b6*(x^4*y^2 + x^2*y^4) ...
   + c6*(x^4 + y^4) ...
   + d6*x^2*y^2 ...
   + e6*(x^2 + y^2);

W6dot = jacobian(W6,[x y])*f;

F6 = [
    sos(W6 - eps_sos*r2^3), ...
    sos(-W6dot - eps_sos*r2^4), ...
    replace(W6,[x y],[1 0]) == 1, ...
    replace(W6,[x y],[1 1]) >= 1
    ];

obj6 = -(a6 + b6);

sol6 = solvesos(F6,obj6,opts,[a6;b6;c6;d6;e6]);

fprintf('Nonhomogeneous sextic status code: %d\n', sol6.problem);
fprintf('Nonhomogeneous sextic solver info: %s\n', sol6.info);

nonhomogeneous_sextic_feasible = (sol6.problem == 0);

if nonhomogeneous_sextic_feasible
    a6sol = value(a6);
    b6sol = value(b6);
    c6sol = value(c6);
    d6sol = value(d6);
    e6sol = value(e6);

    fprintf('Nonhomogeneous sextic SOS search feasible.\n');
    fprintf('Recovered sextic Lyapunov function:\n');
    fprintf(['  W6(x,y) = %.12g*(x^6+y^6) + %.12g*(x^4*y^2+x^2*y^4)' ...
             ' + %.12g*(x^4+y^4) + %.12g*x^2*y^2 + %.12g*(x^2+y^2)\n'], ...
             a6sol,b6sol,c6sol,d6sol,e6sol);

    sextic_mass = abs(a6sol) + abs(b6sol);
    fprintf('Total sextic contribution |a6|+|b6| = %.12g\n\n', sextic_mass);

    if sextic_mass < 1e-8
        warning('The returned nonhomogeneous sextic is effectively lower-order.');
    end
else
    fprintf('Nonhomogeneous sextic SOS search infeasible or numerically rejected.\n\n');
end

%% Analytical check

fprintf('------------------------------------------------------------\n');
fprintf('Analytical check of exact quartic against SOS margins\n');
fprintf('------------------------------------------------------------\n');

fprintf('For V = x^4+y^4, one has Vdot = -4 sin(theta) (x^6+y^6).\n');
fprintf('Also x^4+y^4 >= (1/2)(x^2+y^2)^2 and\n');
fprintf('x^6+y^6 >= (1/4)(x^2+y^2)^3.\n');
fprintf('Thus V - eps*(x^2+y^2)^2 is nonnegative if eps <= 1/2,\n');
fprintf('and -Vdot - eps*(x^2+y^2)^3 is nonnegative if eps <= sin(theta).\n\n');

%% Grids for plotting

xx = linspace(xmin-0.2,xmax+0.2,Nx);
yy = linspace(ymin-0.2,ymax+0.2,Ny);
[X,Y] = meshgrid(xx,yy);

Vexact_grid = X.^4 + Y.^4;
W4_grid = a4sol*(X.^4 + Y.^4) + b4sol*(X.^2).*(Y.^2);

if homogeneous_sextic_feasible
    H6_grid = ah6sol*(X.^6 + Y.^6) ...
        + bh6sol*((X.^4).*(Y.^2) + (X.^2).*(Y.^4));
end

if nonhomogeneous_sextic_feasible
    W6_grid = a6sol*(X.^6 + Y.^6) ...
        + b6sol*((X.^4).*(Y.^2) + (X.^2).*(Y.^4)) ...
        + c6sol*(X.^4 + Y.^4) ...
        + d6sol*(X.^2).*(Y.^2) ...
        + e6sol*(X.^2 + Y.^2);
end

%% Trajectories

ICs = [
    -1.5  1.2;
    -1.3  0.6;
    -1.2 -1.0;
    -0.8  1.4;
     1.5  1.2;
     1.3  0.6;
     1.2 -1.0;
     0.8  1.4;
    -1.4  0.0;
     1.4  0.0;
     0.0  1.4;
     0.0 -1.4];

odeopts = odeset('RelTol',1e-9,'AbsTol',1e-11);
tspan = [0 60];

Traj = cell(size(ICs,1),1);

for k = 1:size(ICs,1)
    [~,sol] = ode45(@(t,z) ode_system(t,z,theta), tspan, ICs(k,:)', odeopts);
    Traj{k} = sol;
end

%% Figure 1: exact quartic surface

fig1 = figure;
set(fig1,'Units','centimeters','Position',[2 2 12 11]);

surf(X,Y,Vexact_grid,'EdgeColor','none');
hold on
contour3(X,Y,Vexact_grid,14,'LineColor',[0.85 0.85 0.85],'LineWidth',0.5);
colormap(parula)
camlight
lighting gouraud

set(gca,'FontSize',24)
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold')
zlabel('$V(x,y)$','Interpreter','latex','FontSize',24,'FontWeight','bold')
box on
grid on
view(42,28)

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.10 0.12 0.82 0.80];

if saveFigures
    save_pdf_figure(fig1,'homogeneous_cubic_exact_quartic_surface.pdf');
end

%% Figure 2: exact quartic contours and trajectories

fig2 = figure;
set(fig2,'Units','centimeters','Position',[2 2 12 11]);

contourf(X,Y,Vexact_grid,16,'LineColor','none');
colormap(parula)
hold on
contour(X,Y,Vexact_grid,12,'LineColor',[0.85 0.85 0.85],'LineWidth',0.4);

nq = 19;
[xq,yq] = meshgrid(linspace(xmin,xmax,nq),linspace(ymin,ymax,nq));
uq = -s*xq.^3 + c*yq.^3;
vq = -c*xq.^3 - s*yq.^3;

quiver(xq,yq,uq,vq,0.8,'Color',[0.55 0.55 0.55],'LineWidth',0.9);

for k = 1:length(Traj)
    sol = Traj{k};
    plot(sol(:,1),sol(:,2),'-','Color',red,'LineWidth',lw_main);
    plot(sol(1,1),sol(1,2),'o', ...
        'Color',black, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',7, ...
        'LineWidth',1.2);
end

plot(0,0,'o', ...
    'Color',green, ...
    'MarkerFaceColor',green, ...
    'MarkerSize',16, ...
    'LineWidth',1.5);

set(gca,'FontSize',24)
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold');
box on
axis equal
axis([xmin xmax ymin ymax])

if saveFigures
    save_pdf_figure(fig2,'quartic_sos_contours_trajectories.pdf');
end

%% Figure 3: exact quartic versus quartic SOS

fig3 = figure;
set(fig3,'Units','centimeters','Position',[2 2 12 11]);

contour(X,Y,Vexact_grid,14,'Color',blue,'LineWidth',lw_main);
hold on
contour(X,Y,W4_grid,14,'LineStyle','--','Color',purple,'LineWidth',lw_main);
contour(X,Y,Vexact_grid,12,'LineColor',[0.85 0.85 0.85],'LineWidth',0.4);

for k = 1:length(Traj)
    sol = Traj{k};
    plot(sol(:,1),sol(:,2),'-','Color',[0.7 0.7 0.7],'LineWidth',1.2);
end

plot(0,0,'o', ...
    'Color',green, ...
    'MarkerFaceColor',green, ...
    'MarkerSize',16, ...
    'LineWidth',1.5);

set(gca,'FontSize',24)
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold');
box on
axis equal
axis([xmin xmax ymin ymax])

if saveFigures
    save_pdf_figure(fig3,'quartic_exact_vs_sos_contours.pdf');
end

%% Figure 4: quartic SOS versus nonhomogeneous sextic SOS

if nonhomogeneous_sextic_feasible
    fig4 = figure;
    set(fig4,'Units','centimeters','Position',[2 2 12 11]);

    contour(X,Y,W4_grid,14,'Color',blue,'LineWidth',lw_main);
    hold on
    contour(X,Y,W6_grid,14,'LineStyle','--','Color',purple,'LineWidth',lw_main);

    for k = 1:length(Traj)
        sol = Traj{k};
        plot(sol(:,1),sol(:,2),'-','Color',[0.7 0.7 0.7],'LineWidth',1.0);
    end

    plot(0,0,'o', ...
        'Color',green, ...
        'MarkerFaceColor',green, ...
        'MarkerSize',10, ...
        'LineWidth',1.5);

    set(gca,'FontSize',24)
    xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
    ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold')
    axis equal
    axis([xmin xmax ymin ymax])
    grid on
    box on

    ax = gca;
    ax.Units = 'normalized';
    ax.Position = [0.10 0.12 0.82 0.80];

    if saveFigures
        save_pdf_figure(fig4,'homogeneous_cubic_quartic_vs_sextic_sos.pdf');
    end
end

%% Figure 5: compare quartic, homogeneous sextic, and nonhomogeneous sextic

if homogeneous_sextic_feasible && nonhomogeneous_sextic_feasible
    fig5 = figure;
    set(fig5,'Units','centimeters','Position',[2 2 12 11]);

    contour(X,Y,W4_grid,12,'Color',blue,'LineWidth',lw_main);
    hold on
    contour(X,Y,H6_grid,12,'LineStyle','-.','Color',red,'LineWidth',lw_main);
    contour(X,Y,W6_grid,12,'LineStyle','--','Color',purple,'LineWidth',lw_main);

    set(gca,'FontSize',24)
    xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold')
    ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold')
    axis equal
    axis([xmin xmax ymin ymax])
    grid on
    box on

    ax = gca;
    ax.Units = 'normalized';
    ax.Position = [0.10 0.12 0.82 0.80];

    if saveFigures
        save_pdf_figure(fig5,'homogeneous_cubic_all_sos_contours.pdf');
    end
end

%% Figure 6: decay of the exact quartic along one trajectory

[tcheck,solcheck] = ode45(@(t,z) ode_system(t,z,theta), ...
    [0 60], [1.4;1.0], odeopts);

Vtraj_exact = solcheck(:,1).^4 + solcheck(:,2).^4;

fig6 = figure;
set(fig6,'Units','centimeters','Position',[2 2 12 9]);

plot(tcheck,Vtraj_exact,'Color',blue,'LineWidth',lw_main);

set(gca,'FontSize',24)
xlabel('$t$','Interpreter','latex','FontSize',24,'FontWeight','bold')
ylabel('$V_{\mathrm{exact}}(x(t),y(t))$', ...
    'Interpreter','latex','FontSize',24,'FontWeight','bold')
grid on
box on

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.12 0.16 0.82 0.76];

if saveFigures
    save_pdf_figure(fig6,'homogeneous_cubic_exact_quartic_decay.pdf');
end

fprintf('Finished. Figures saved to current directory.\n');

%% Local functions

function dz = ode_system(~,z,theta)

    x = z(1);
    y = z(2);

    s = sin(theta);
    c = cos(theta);

    dz = [-s*x^3 + c*y^3;
          -c*x^3 - s*y^3];
end

function save_pdf_figure(fig,fileName)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n', fileName);
end