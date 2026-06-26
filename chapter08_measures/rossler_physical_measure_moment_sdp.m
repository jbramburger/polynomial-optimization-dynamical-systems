%% rossler_physical_measure_moment_sdp.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Invariant Measures, Ergodic Optimization, and Duality
% Example: Moment-SDP approximation of a physical measure for the Rossler system
%
% This script approximates the physical invariant measure of the Rossler
% system
%
%     x' = -y - z,
%     y' =  x + a y,
%     z' =  b + z(x-c),
%
% with parameters a = b = 0.1 and c = 18.
%
% The script:
%
%   1. Simulates a long trajectory on the chaotic attractor.
%   2. Computes empirical first moments from the trajectory.
%   3. Solves a moment SDP on a rescaled unit ball.
%   4. Fits an invariant measure whose first moments match the data.
%   5. Reconstructs a polynomial density from the moment sequence.
%   6. Plots the attractor and component time series colored by density.
%
% Pressing Run should reproduce the numerical output and save:
%
%     rossler_density_attractor.png
%     rossler_density_components.pdf
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

a = 0.1;
b = 0.1;
c = 18;

momentDegree = 10;

Ttrans = 300;
Tdata  = 2000;
dt     = 0.02;

x0 = [1; 1; 1];

ballPadding = 1.10;

plotStride = 1;
componentPlotTime = 300;

solver_name = 'mosek';
verbose_solver = 1;

saveFigures = true;

rng(1);

%% Solver options

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver);

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.black  = [0, 0, 0];

%% Simulate Rossler trajectory

fprintf('\n============================================================\n');
fprintf('Rossler physical measure via a moment SDP\n');
fprintf('============================================================\n\n');

fprintf('Parameters:\n');
fprintf('  a = %.6g\n', a);
fprintf('  b = %.6g\n', b);
fprintf('  c = %.6g\n', c);
fprintf('Moment degree d = %d\n\n', momentDegree);

fprintf('Simulating transient...\n');

odeopts = odeset('RelTol',1e-9,'AbsTol',1e-11);

[~, Xtrans] = ode45(@(t,x) rossler_rhs(t,x,a,b,c), ...
    [0 Ttrans], x0, odeopts);

xstart = Xtrans(end,:).';

fprintf('Simulating data trajectory...\n');

tspan = 0:dt:Tdata;

[t, X] = ode45(@(t,x) rossler_rhs(t,x,a,b,c), ...
    tspan, xstart, odeopts);

xdat = X(:,1);
ydat = X(:,2);
zdat = X(:,3);

%% Affine rescaling to the unit ball

center = mean(X,1).';

distances = vecnorm(X - center.', 2, 2);
R = ballPadding*max(distances);

Sdata = (X - center.')/R;

Xsc = Sdata(:,1);
Ysc = Sdata(:,2);
Zsc = Sdata(:,3);

fprintf('\nRescaling data to unit ball:\n');
fprintf('  center = [%.8f %.8f %.8f]\n', center(1), center(2), center(3));
fprintf('  R      = %.8f\n', R);
fprintf('  max scaled radius = %.8f\n', max(vecnorm(Sdata,2,2)));

%% Empirical first moments in scaled variables

empiricalFirstMoments = mean(Sdata,1).';

fprintf('\nEmpirical first moments in scaled variables:\n');
fprintf('  E[X] = %.8f\n', empiricalFirstMoments(1));
fprintf('  E[Y] = %.8f\n', empiricalFirstMoments(2));
fprintf('  E[Z] = %.8f\n', empiricalFirstMoments(3));

%% Solve moment-fitting SDP

fprintf('\nSolving degree-%d moment SDP...\n', momentDegree);

[yval, basis] = solve_rossler_physical_moment_sdp( ...
    momentDegree, a, b, c, center, R, empiricalFirstMoments, opts);

%% Reconstruct polynomial density relative to Lebesgue measure

fprintf('\nReconstructing polynomial density on the unit ball...\n');

[densityCoeff, densityBasis] = reconstruct_density_3d( ...
    yval, basis, momentDegree);

%% Evaluate reconstructed density along trajectory

rho = zeros(size(Xsc));

for k = 1:length(rho)
    rho(k) = eval_poly_coeff_3d( ...
        densityCoeff, densityBasis, Xsc(k), Ysc(k), Zsc(k));
end

% Polynomial density reconstructions may have small oscillatory negative
% values. For visualization, clip at zero and color on a logarithmic scale.
rhoClip = max(rho,0);

positiveVals = rhoClip(rhoClip > 0);

if isempty(positiveVals)
    rhoFloor = 1e-12;
else
    rhoFloor = max(prctile(positiveVals,1),1e-12);
end

rhoColor = log10(max(rhoClip,rhoFloor));

fprintf('\nDensity summary along trajectory:\n');
fprintf('  min rho         = %.4e\n', min(rho));
fprintf('  max rho         = %.4e\n', max(rho));
fprintf('  min clipped rho = %.4e\n', min(rhoClip));
fprintf('  max clipped rho = %.4e\n', max(rhoClip));

%% Subsample for plotting

idxp = 1:plotStride:length(t);

tp = t(idxp);
xp = xdat(idxp);
yp = ydat(idxp);
zp = zdat(idxp);
cp = rhoColor(idxp);

%% Figure 1: attractor colored by reconstructed density

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 13 11]);

surface( ...
    [xp xp], [yp yp], [zp zp], [cp cp], ...
    'FaceColor','none', ...
    'EdgeColor','interp', ...
    'LineWidth',1.5);

colormap(parula);
colorbar;

view(-40,25);
axis equal tight;
grid on
box on

set(gca,'FontSize',42,'TickLabelInterpreter','latex');

xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');
zlabel('$z$','Interpreter','latex','FontSize',42,'FontWeight','bold');

if saveFigures
    exportgraphics(fig1,'rossler_density_attractor.png', ...
        'Resolution',600, ...
        'BackgroundColor','white');
    fprintf('Saved figure: rossler_density_attractor.png\n');
end

%% Figure 2: components colored by reconstructed density

keep = tp <= componentPlotTime;

tt = tp(keep);
xx = xp(keep);
yy = yp(keep);
zz = zp(keep);
cc = cp(keep);

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 13 11]);

tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
surface([tt tt],[xx xx],zeros(length(tt),2),[cc cc], ...
    'FaceColor','none', ...
    'EdgeColor','interp', ...
    'LineWidth',5);
view(2);
grid on
box on
ylabel('$x(t)$','Interpreter','latex','FontSize',32,'FontWeight','bold');
set(gca,'FontSize',32,'TickLabelInterpreter','latex');
colorbar;

nexttile;
surface([tt tt],[yy yy],zeros(length(tt),2),[cc cc], ...
    'FaceColor','none', ...
    'EdgeColor','interp', ...
    'LineWidth',5);
view(2);
grid on
box on
ylabel('$y(t)$','Interpreter','latex','FontSize',32,'FontWeight','bold');
set(gca,'FontSize',32,'TickLabelInterpreter','latex');
colorbar;

nexttile;
surface([tt tt],[zz zz],zeros(length(tt),2),[cc cc], ...
    'FaceColor','none', ...
    'EdgeColor','interp', ...
    'LineWidth',5);
view(2);
grid on
box on
xlabel('$t$','Interpreter','latex','FontSize',32,'FontWeight','bold');
ylabel('$z(t)$','Interpreter','latex','FontSize',32,'FontWeight','bold');
set(gca,'FontSize',32,'TickLabelInterpreter','latex');
colorbar;

if saveFigures
    save_pdf_figure(fig2,'rossler_density_components.pdf');
end

fprintf('\nFinished Rossler physical-measure moment-SDP example.\n');

%% Local functions

function dx = rossler_rhs(~,x,a,b,c)

    dx = zeros(3,1);

    dx(1) = -x(2) - x(3);
    dx(2) =  x(1) + a*x(2);
    dx(3) =  b + x(3)*(x(1) - c);
end

function [yval, basis] = solve_rossler_physical_moment_sdp( ...
    d, a, b, c, center, R, targetFirstMoments, opts)

    basis = monomial_basis_3d(d);
    N = size(basis,1);

    y = sdpvar(N,1);

    idx = make_index_map_3d(basis);

    constraints = [];

    %% Probability normalization

    constraints = [constraints, y(get_idx_3d(idx,0,0,0)) == 1];

    %% Positivity of measure: moment matrix

    rMoment = floor(d/2);
    basisMoment = monomial_basis_3d(rMoment);

    M = build_moment_matrix_3d(y,idx,basisMoment,[0 0 0]);
    constraints = [constraints, M >= 0];

    %% Support on the unit ball

    rLocalizing = floor((d - 2)/2);

    if rLocalizing >= 0

        basisLocalizing = monomial_basis_3d(rLocalizing);

        Mg = build_ball_localizing_matrix_3d(y,idx,basisLocalizing);
        constraints = [constraints, Mg >= 0];
    end

    %% Scaled Rossler vector field in variables X,Y,Z
    %
    % Original variables:
    %
    %     q = center + R s,
    %
    % where s = (X,Y,Z). The scaled dynamics are
    %
    %     s' = f(center + R s)/R.

    cx = center(1);
    cy = center(2);
    cz = center(3);

    % F1 = (-y-z)/R = const - Y - Z.
    F1_terms = [
        -(cy + cz)/R, 0, 0, 0;
        -1,           0, 1, 0;
        -1,           0, 0, 1
    ];

    % F2 = (x + a y)/R = const + X + aY.
    F2_terms = [
        (cx + a*cy)/R, 0, 0, 0;
        1,             1, 0, 0;
        a,             0, 1, 0
    ];

    % F3 = (b + z(x-c))/R
    %    = const + cz X + (cx-c) Z + R XZ.
    F3_terms = [
        (b + cz*(cx-c))/R, 0, 0, 0;
        cz,                1, 0, 0;
        (cx-c),            0, 0, 1;
        R,                 1, 0, 1
    ];

    %% Weak invariance constraints
    %
    %     int grad(p) . F dmu = 0
    %
    % for all test polynomials p with degree at most d - 1.

    degf = 2;
    testDegree = d + 1 - degf;

    if testDegree >= 1

        basisTest = monomial_basis_3d(testDegree);

        for k = 1:size(basisTest,1)

            alpha = basisTest(k,:);

            ax = alpha(1);
            ay = alpha(2);
            az = alpha(3);

            if ax == 0 && ay == 0 && az == 0
                continue;
            end

            expr = 0;

            if ax > 0
                expr = expr + ax*apply_terms_to_moment( ...
                    y,idx,ax-1,ay,az,F1_terms);
            end

            if ay > 0
                expr = expr + ay*apply_terms_to_moment( ...
                    y,idx,ax,ay-1,az,F2_terms);
            end

            if az > 0
                expr = expr + az*apply_terms_to_moment( ...
                    y,idx,ax,ay,az-1,F3_terms);
            end

            constraints = [constraints, expr == 0]; %#ok<AGROW>
        end
    end

    %% Fit empirical first moments

    EX = moment_expr_3d(y,idx,1,0,0);
    EY = moment_expr_3d(y,idx,0,1,0);
    EZ = moment_expr_3d(y,idx,0,0,1);

    objective = (EX - targetFirstMoments(1))^2 ...
              + (EY - targetFirstMoments(2))^2 ...
              + (EZ - targetFirstMoments(3))^2;

    diagnostics = optimize(constraints, objective, opts);

    if diagnostics.problem ~= 0
        warning('Moment-fitting SDP diagnostic: %s', diagnostics.info);
    end

    fprintf('moment-fitting objective = %.4e\n', value(objective));
    fprintf('matched E[X] = %.8f, target %.8f\n', value(EX), targetFirstMoments(1));
    fprintf('matched E[Y] = %.8f, target %.8f\n', value(EY), targetFirstMoments(2));
    fprintf('matched E[Z] = %.8f, target %.8f\n', value(EZ), targetFirstMoments(3));

    yval = value(y);
end

function expr = apply_terms_to_moment(y, idx, a0, b0, c0, terms)

    expr = 0;

    for r = 1:size(terms,1)

        coeff = terms(r,1);
        da = terms(r,2);
        db = terms(r,3);
        dc = terms(r,4);

        expr = expr + coeff*moment_expr_3d( ...
            y,idx,a0+da,b0+db,c0+dc);
    end
end

function basis = monomial_basis_3d(d)

    basis = [];

    for total = 0:d
        for ax = 0:total
            for ay = 0:(total-ax)

                az = total - ax - ay;
                basis = [basis; ax ay az]; %#ok<AGROW>
            end
        end
    end
end

function idx = make_index_map_3d(basis)

    idx = containers.Map();

    for k = 1:size(basis,1)

        key = sprintf('%d_%d_%d', ...
            basis(k,1), basis(k,2), basis(k,3));

        idx(key) = k;
    end
end

function k = get_idx_3d(idx, a, b, c)

    key = sprintf('%d_%d_%d', a, b, c);

    if ~isKey(idx,key)
        error('Moment y_(%d,%d,%d) is outside the truncation.', a, b, c);
    end

    k = idx(key);
end

function expr = moment_expr_3d(y, idx, a, b, c)

    expr = y(get_idx_3d(idx,a,b,c));
end

function M = build_moment_matrix_3d(y, idx, basisMoment, shift)

    n = size(basisMoment,1);

    M = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            aa = basisMoment(i,1) + basisMoment(j,1) + shift(1);
            bb = basisMoment(i,2) + basisMoment(j,2) + shift(2);
            cc = basisMoment(i,3) + basisMoment(j,3) + shift(3);

            M(i,j) = moment_expr_3d(y,idx,aa,bb,cc);
            M(j,i) = M(i,j);
        end
    end
end

function Mg = build_ball_localizing_matrix_3d(y, idx, basisLocalizing)

    n = size(basisLocalizing,1);

    Mg = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            aa = basisLocalizing(i,1) + basisLocalizing(j,1);
            bb = basisLocalizing(i,2) + basisLocalizing(j,2);
            cc = basisLocalizing(i,3) + basisLocalizing(j,3);

            val = moment_expr_3d(y,idx,aa,bb,cc) ...
                - moment_expr_3d(y,idx,aa+2,bb,cc) ...
                - moment_expr_3d(y,idx,aa,bb+2,cc) ...
                - moment_expr_3d(y,idx,aa,bb,cc+2);

            Mg(i,j) = val;
            Mg(j,i) = val;
        end
    end
end

function [coeff, basis] = reconstruct_density_3d(yval, basisY, d)
% Reconstruct rho_d(X,Y,Z) = coeff' z_d(X,Y,Z) satisfying
%
%     int_B z_d rho_d dX dY dZ = y,
%
% where B is the unit ball in R^3.

    basis = monomial_basis_3d(d);
    N = size(basis,1);

    massMatrix = zeros(N,N);
    rhs = zeros(N,1);

    idxY = make_index_map_3d(basisY);

    for i = 1:N

        aa = basis(i,1);
        bb = basis(i,2);
        cc = basis(i,3);

        rhs(i) = yval(get_idx_3d(idxY,aa,bb,cc));

        for j = 1:N
            pp = basis(i,:) + basis(j,:);
            massMatrix(i,j) = ball_moment_3d(pp(1),pp(2),pp(3));
        end
    end

    reg = 1e-12;
    coeff = (massMatrix + reg*eye(N))\rhs;
end

function val = ball_moment_3d(a,b,c)
% Integral over the unit ball in R^3 of X^a Y^b Z^c dX dY dZ.

    if mod(a,2) == 1 || mod(b,2) == 1 || mod(c,2) == 1
        val = 0;
        return;
    end

    i = a/2;
    j = b/2;
    k = c/2;

    val = gamma(i+0.5)*gamma(j+0.5)*gamma(k+0.5) ...
        / gamma(i+j+k+2.5);
end

function val = eval_poly_coeff_3d(coeff, basis, x, y, z)

    val = 0;

    for k = 1:length(coeff)

        val = val ...
            + coeff(k)*x^basis(k,1)*y^basis(k,2)*z^basis(k,3);
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