%% motzkin.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Polynomial Optimization and Sum-of-Squares Programming
% Example: Global SOS versus Putinar certificate for the Motzkin polynomial
%
% This script compares three SOS-based computations for the Motzkin
% polynomial
%
%     m(x,y) = x^4 y^2 + x^2 y^4 + 1 - 3 x^2 y^2.
%
% The script:
%
%   1. Tests whether m(x,y) is globally SOS.
%   2. Attempts the global SOS lower-bound problem and verifies that it is
%      infeasible: m(x,y) - gamma is not SOS for any constant gamma.
%   3. Computes a Putinar-type SOS lower bound on the unit disk.
%   4. Produces the contour figure used in the book.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     motzkin.pdf
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
figureFile = 'motzkin.pdf';

opts = sdpsettings('solver','mosek','verbose',0);

%% Book plotting style

blue   = [0, 92, 175]/255;
red    = [200, 50, 50]/255;
green  = [0, 140, 90]/255;

%% Define the Motzkin polynomial

m_fun = @(x,y) x.^4 .* y.^2 + x.^2 .* y.^4 + 1 - 3*x.^2 .* y.^2;

sdpvar x y gamma
m = x^4*y^2 + x^2*y^4 + 1 - 3*x^2*y^2;

%% 1. Direct global SOS feasibility test

sol_zero = solvesos(sos(m), [], opts);

%% 2. Attempt a global SOS lower bound

% No artificial lower bound is imposed on gamma. In particular, a value such
% as gamma = -10 must not be mistaken for a feasible SOS lower bound after an
% unsuccessful solver call.
Fglob = sos(m - gamma);
sol_glob = solvesos(Fglob, -gamma, opts, gamma);

global_sos_feasible = (sol_glob.problem == 0);

if global_sos_feasible
    gamma_global = value(gamma);
else
    gamma_global = NaN;
end

%% 3. Putinar lower bound on the unit disk

g = 1 - x^2 - y^2;

[s0,c0] = polynomial([x y], 6);
[s1,c1] = polynomial([x y], 4);

expr = m - gamma - s0 - s1*g;

Fdisk = [
    sos(s0), ...
    sos(s1), ...
    coefficients(expr,[x y]) == 0, ...
    -10 <= gamma <= 1
    ];

sol_disk = solvesos(Fdisk, -gamma, opts, [gamma; c0; c1]);

disk_sos_feasible = (sol_disk.problem == 0);

if disk_sos_feasible
    gamma_disk = value(gamma);
else
    gamma_disk = NaN;
    warning('Putinar optimization failed: %s', sol_disk.info);
end

%% Exact values and minimizers

true_global_min = 0;
true_disk_min   = 0.5;

global_min_pts = [ 1,  1;
                  -1,  1;
                   1, -1;
                  -1, -1];

c = 1/sqrt(2);
disk_min_pts = [ c,  c;
                -c,  c;
                 c, -c;
                -c, -c];

%% Print results

fprintf('\n============================================================\n');
fprintf('Motzkin polynomial: global SOS versus Putinar on the disk\n');
fprintf('============================================================\n\n');

fprintf('True global minimum on R^2:        %.12f\n', true_global_min);
if global_sos_feasible
    fprintf('Computed global SOS lower bound:  %.12f\n', gamma_global);
    fprintf('Gap:                              %.12e\n\n', ...
        true_global_min - gamma_global);
else
    fprintf('Global SOS lower-bound problem:    infeasible\n');
    fprintf('No constant gamma makes m-gamma SOS.\n');
    fprintf('  solver status code: %d\n', sol_glob.problem);
    fprintf('  solver message:     %s\n\n', sol_glob.info);
end

fprintf('Direct SOS feasibility test for m(x,y):\n');
fprintf('  solver status code: %d\n', sol_zero.problem);
fprintf('  solver message:     %s\n\n', sol_zero.info);

fprintf('True minimum on the unit disk:     %.12f\n', true_disk_min);
if disk_sos_feasible
    fprintf('Computed Putinar lower bound:      %.12f\n', gamma_disk);
    fprintf('Gap:                               %.12e\n\n', ...
        true_disk_min - gamma_disk);
else
    fprintf('Putinar lower-bound problem:       solver failure\n');
    fprintf('No numerical bound is reported.\n\n');
end

fprintf('Global minimizers:\n');
disp(global_min_pts);

fprintf('Constrained minimizers on the unit disk:\n');
disp(disk_min_pts);

%% Produce contour figure

xx = linspace(-1.5, 1.5, 500);
yy = linspace(-1.5, 1.5, 500);
[X,Y] = meshgrid(xx,yy);
M = m_fun(X,Y);

fig = figure;
set(fig,'Units','centimeters','Position',[2 2 12 11]);

contourf(X,Y,M,16,'LineColor','none');
colormap(parula)
caxis([min(M(:)), max(M(:))*0.2])
hold on

contour(X,Y,M,10,'LineColor',[0.85 0.85 0.85],'LineWidth',0.4);

th = linspace(0,2*pi,400);
plot(cos(th), sin(th), '--', ...
    'Color', red, ...
    'LineWidth', 3);

plot(global_min_pts(:,1), global_min_pts(:,2), ...
    'o', ...
    'Color', blue, ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 16, ...
    'LineWidth', 2);

plot(disk_min_pts(:,1), disk_min_pts(:,2), ...
    's', ...
    'Color', green, ...
    'MarkerFaceColor', green, ...
    'MarkerSize', 16, ...
    'LineWidth', 1.8);

set(gca, 'FontSize', 24)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')

axis equal
axis([-1.5 1.5 -1.5 1.5])
grid on
box on
colorbar

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.08 0.22 0.84 0.74];

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
