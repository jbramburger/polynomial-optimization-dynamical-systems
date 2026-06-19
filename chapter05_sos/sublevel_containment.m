%% sublevel_containment.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Polynomial Optimization and Sum-of-Squares Programming
% Example: Sublevel set containment using SOS certificates
%
% This script certifies the containment
%
%     { (x,y) : p(x,y) <= c } subset { (x,y) : q(x,y) <= 1 }
%
% inside the ambient ball
%
%     K = { (x,y) : x^2 + y^2 <= 2 },
%
% where
%
%     p(x,y) = x^2 + y^2.
%
% The largest certified value of c is found by bisection. For each fixed
% c, feasibility is checked using a Putinar-type SOS certificate.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     sublevel_containment.pdf
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

saveFigure = true;
figureFile = 'sublevel_containment.pdf';

opts = sdpsettings('solver','mosek','verbose',0);

%% Book plotting style

blue   = [0, 92, 175]/255;
red    = [200, 50, 50]/255;

%% Polynomial definitions

p_fun = @(x,y) x.^2 + y.^2;

q_fun = @(x,y) 0.6*x.^4 + y.^4 + 0.5*x.^2.*y.^2 + 0.25*x.^3 ...
             - 0.2*x.*y + 0.35*x.^2 + 0.8*y.^2;

sdpvar x y

p = x^2 + y^2;

q = 0.6*x^4 + y^4 + 0.5*x^2*y^2 + 0.25*x^3 ...
  - 0.2*x*y + 0.35*x^2 + 0.8*y^2;

% Ambient region K = {x^2 + y^2 <= 2}
gK = 2 - x^2 - y^2;

%% Feasibility checker for fixed c

is_feasible = @(cval) check_feasibility(cval, x, y, p, q, gK, opts);

%% Bisection on c

clow  = 0.0;
chigh = 2.0;

% Increase upper endpoint until it is infeasible.
for j = 1:20
    if ~is_feasible(chigh)
        break
    end
    chigh = chigh + 1.0;
end

tol = 1e-6;
maxit = 40;

fprintf('\nBisection search for certified sublevel value c\n');
fprintf('------------------------------------------------------------\n');

for k = 1:maxit
    cmid = 0.5*(clow + chigh);
    feas = is_feasible(cmid);

    fprintf('iter = %2d, c = %.8f, feasible = %d\n', k, cmid, feas);

    if feas
        clow = cmid;
    else
        chigh = cmid;
    end

    if chigh - clow < tol
        break
    end
end

c_sos = clow;

%% Print results

fprintf('\n============================================================\n');
fprintf('Sublevel set containment example\n');
fprintf('============================================================\n\n');

fprintf('Certified value of c:              %.8f\n', c_sos);
fprintf('Radius of certified disk sqrt(c):  %.8f\n\n', sqrt(c_sos));

%% Produce figure

xx = linspace(-1.5,1.5,500);
yy = linspace(-1.5,1.5,500);
[X,Y] = meshgrid(xx,yy);

Q = q_fun(X,Y);

th = linspace(0,2*pi,500);
r  = sqrt(c_sos);
R  = sqrt(2);

fig = figure;
set(fig,'Units','centimeters','Position',[2 2 12 11]);

% Filled contours of q
contourf(X,Y,Q,12,'LineColor','none');
colormap(parula)
caxis([min(Q(:)), max(Q(:))*0.4])
alpha(0.85)
colorbar
hold on

% Certified disk
fill(r*cos(th), r*sin(th), ...
    blue, ...
    'FaceAlpha', 0.08, ...
    'EdgeColor','none');

% Light contour lines
contour(X,Y,Q,6,'LineColor',[0.85 0.85 0.85],'LineWidth',0.5);

% Boundary q = 1
contour(X,Y,Q,[1 1], ...
    'Color', red, ...
    'LineWidth', 3);

% Certified p = c circle
plot(r*cos(th), r*sin(th), ...
    '.', ...
    'Color', 'w', ...
    'LineWidth', 4.5);

% Ambient region boundary
plot(R*cos(th), R*sin(th), '--', ...
    'Color', 'k', ...
    'LineWidth', 2.0);

set(gca, 'FontSize', 18)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')

axis equal
axis([-1.5 1.5 -1.5 1.5])
grid on
box on

%% Save figure

if saveFigure
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'Units', 'centimeters');

    pos = get(fig, 'Position');

    set(fig, 'PaperSize', [pos(3) pos(4)]);
    set(fig, 'PaperPositionMode', 'manual');
    set(fig, 'PaperPosition', [0 0 pos(3) pos(4)]);

    print(fig, figureFile, '-dpdf');
    fprintf('Saved figure: %s\n', figureFile);
end

%% Local function

function feas = check_feasibility(cval, x, y, p, q, gK, opts)

    % SOS certificate:
    %
    %     1 - q = s0 + s1 gK + s2 (c - p),
    %
    % with s0, s1, s2 SOS. This certifies q <= 1 on
    % {gK >= 0, c - p >= 0}.

    [s0,c0] = polynomial([x y],4);
    [s1,c1] = polynomial([x y],2);
    [s2,c2] = polynomial([x y],2);

    expr = 1 - q - s0 - s1*gK - s2*(cval - p);

    F = [
        sos(s0), ...
        sos(s1), ...
        sos(s2), ...
        coefficients(expr,[x y]) == 0
        ];

    sol = solvesos(F, [], opts, [c0; c1; c2]);

    feas = (sol.problem == 0);
end