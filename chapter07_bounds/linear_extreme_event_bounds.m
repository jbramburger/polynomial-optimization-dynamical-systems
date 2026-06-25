%% linear_extreme_event_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Optimal Bounds on Dynamical Systems via Auxiliary Functions
% Example: Finite-time extreme-event bounds for a nonnormal linear system
%
% This script computes SOS upper bounds on the finite-time maximum of x^2
% for the nonnormal linear system
%
%     x' = -x + a y,
%     y' = -y,
%
% with initial conditions in the unit disk.
%
% Pressing Run should reproduce the numerical output and save:
%
%     linear_extreme_event_bounds.csv
%     extreme_phase_plane.pdf
%     extreme_defect_1.pdf
%     extreme_defect_2.pdf
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

a = 8;

Tvals = [0.25 0.5 0.75 1 1.5 2];
degVals = [2 4 6];

verbose_solver = 0;

saveFigures = true;
saveTable = true;

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.orange = [230, 120, 20]/255;
S.purple = [170, 90, 160]/255;
S.cyan   = [0, 170, 200]/255;
S.black  = [0, 0, 0];
S.gray   = [0.40, 0.40, 0.40];
S.fill_blue = [120, 170, 220]/255;

%% Solve SOS bounds

bounds = nan(length(Tvals), length(degVals));
exactVals = nan(length(Tvals),1);

fprintf('\n============================================================\n');
fprintf('Finite-time extreme-event bounds for a nonnormal linear system\n');
fprintf('============================================================\n\n');

fprintf('a = %.6g\n\n', a);

for iT = 1:length(Tvals)

    T = Tvals(iT);
    exactVals(iT) = exact_extreme(a,T);

    R = 1.25*sqrt(1 + a^2*T^2);

    for id = 1:length(degVals)

        degV = degVals(id);

        fprintf('Solving T = %.3f, degree = %d\n', T, degV);

        [Cval, solinfo] = solve_extreme_sos(a,T,degV,R,verbose_solver);

        bounds(iT,id) = Cval;

        if solinfo.problem ~= 0
            warning('Solver issue for T=%.3f, degree=%d: %s', ...
                T, degV, solinfo.info);
        end
    end
end

%% Results table

Results = table(Tvals(:), exactVals(:), ...
    bounds(:,1), bounds(:,2), bounds(:,3), ...
    bounds(:,1)-exactVals(:), ...
    bounds(:,2)-exactVals(:), ...
    bounds(:,3)-exactVals(:), ...
    'VariableNames', {'T','Exact','d2','d4','d6', ...
    'gap_d2','gap_d4','gap_d6'});

disp(Results);

if saveTable
    writetable(Results,'linear_extreme_event_bounds.csv');
    fprintf('Saved table: linear_extreme_event_bounds.csv\n');
end

%% Figure 1: phase-plane trajectories

Tplot = max(Tvals);
tt = linspace(0,Tplot,700);

xmin = -3;
xmax = 3.0;
ymin = -1.2;
ymax = 1.2;

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 12 9]);

hold on
box on

% Background vector field.
nx = 23;
ny = 19;

[Xg,Yg] = meshgrid(linspace(xmin,xmax,nx), linspace(ymin,ymax,ny));

Ug = -Xg + a*Yg;
Vg = -Yg;

Spd = sqrt(Ug.^2 + Vg.^2);
Ug = Ug ./ max(Spd,1e-12);
Vg = Vg ./ max(Spd,1e-12);

quiver(Xg,Yg,Ug,Vg,0.32, ...
    'Color',0.82*[1 1 1], ...
    'LineWidth',0.65, ...
    'MaxHeadSize',0.55);

% Initial set.
theta = linspace(0,2*pi,600);

fill(cos(theta),sin(theta),S.fill_blue, ...
    'FaceAlpha',0.16, ...
    'EdgeColor','none');

h_X0 = plot(cos(theta),sin(theta),'--', ...
    'Color',S.blue, ...
    'LineWidth',2.0);

% Sample trajectories.
ths = linspace(0,2*pi,25);
ths(end) = [];

traj_cols = [S.orange; S.purple; S.green; S.cyan; S.blue];

for k = 1:length(ths)

    th = ths(k);

    x0 = cos(th);
    y0 = sin(th);

    xtraj = exp(-tt).*(x0 + a*tt*y0);
    ytraj = exp(-tt).*y0;

    col = interp1( ...
        linspace(1,length(ths),size(traj_cols,1)), ...
        traj_cols,k,'linear');

    plot(xtraj,ytraj, ...
        'Color',0.75*col + 0.25*S.gray, ...
        'LineWidth',1.25);

    idx = round(0.82*length(tt));
    add_arrow_head_on_curve(xtraj,ytraj,idx,0.045,0.75*col + 0.25*S.gray);
end

h_ic = scatter(cos(ths),sin(ths),46,S.blue,'filled', ...
    'MarkerEdgeColor',0.65*S.blue, ...
    'LineWidth',0.6);

% Extreme trajectory.
tstar = maximizing_time(a,Tplot);

x0s = 1/sqrt(1+a^2*tstar^2);
y0s = a*tstar/sqrt(1+a^2*tstar^2);

xstar = exp(-tt).*(x0s + a*tt*y0s);
ystar = exp(-tt).*y0s;

h_ext = plot(xstar,ystar, ...
    'Color',S.red, ...
    'LineWidth',3.0);

arrow_ids = round([0.35 0.62 0.88]*length(tt));

for jj = arrow_ids
    add_arrow_head_on_curve(xstar,ystar,jj,0.065,S.red);
end

h_maxic = scatter(x0s,y0s,95,S.red,'filled', ...
    'MarkerEdgeColor',0.55*S.red, ...
    'LineWidth',0.8);

h_origin = scatter(0,0,100,S.black,'filled');

axis([xmin xmax ymin ymax]);

xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$y$','Interpreter','latex','FontSize',24,'FontWeight','bold');

grid on

set(gca, ...
    'FontSize',24, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.0, ...
    'GridAlpha',0.16, ...
    'MinorGridAlpha',0.10);

h_traj = plot(nan,nan,'-', ...
    'Color',S.gray, ...
    'LineWidth',1.4);

legend([h_X0,h_ic,h_traj,h_ext,h_maxic,h_origin], ...
    {'$X_0$ (unit disk)', ...
     'sample initial conditions', ...
     'sample trajectories', ...
     'extreme trajectory', ...
     'maximizing initial condition', ...
     'origin'}, ...
     'Interpreter','latex', ...
     'Location','southeast', ...
     'FontSize',18, ...
     'Box','on');

if saveFigures
    save_pdf_figure(fig1,'extreme_phase_plane.pdf');
end

%% Defect visualization

Tvis = max(Tvals);
degVis = max(degVals);
Rvis = 1.25*sqrt(1 + a^2*Tvis^2);

fprintf('\nRe-solving for defect plots: T = %.3f, degree = %d\n', Tvis, degVis);

[~, solinfo, aux] = solve_extreme_sos(a,Tvis,degVis,Rvis,verbose_solver);

if solinfo.problem ~= 0
    warning('Solver issue in defect solve: %s', solinfo.info);
end

tphys = maximizing_time(a,Tvis);
svis = tphys/Tvis;

xgrid = linspace(xmin,xmax,180);
ygrid = linspace(ymin,ymax,180);
[X,Y] = meshgrid(xgrid,ygrid);

D1vals = nan(size(X));
D2vals = nan(size(X));

for k = 1:numel(X)

    if X(k)^2 + Y(k)^2 <= Rvis^2
        D1vals(k) = eval_poly(aux.D1num, aux.vars, [X(k),Y(k),svis]);
        D2vals(k) = eval_poly(aux.D2num, aux.vars, [X(k),Y(k),svis]);
    end
end

D1vals = max(D1vals,1e-12);
D2vals = max(D2vals,1e-12);

ttDef = linspace(0,Tvis,500);

xstarDef = exp(-ttDef).*(x0s + a*ttDef*y0s);
ystarDef = exp(-ttDef).*y0s;

[~,idxStar] = min(abs(ttDef-tphys));

plot_defect_field(X,Y,D1vals,theta,xstarDef,ystarDef, ...
    idxStar,x0s,y0s,S,'extreme_defect_1.pdf',saveFigures);

plot_defect_field(X,Y,D2vals,theta,xstarDef,ystarDef, ...
    idxStar,x0s,y0s,S,'extreme_defect_2.pdf',saveFigures);

fprintf('\nFinished linear extreme-event example.\n');

%% Local functions

function [Cval, sol, aux] = solve_extreme_sos(a,T,degV,R,verbose)

    yalmip('clear');

    sdpvar x y s C

    f1 = -x + a*y;
    f2 = -y;
    Phi = x^2;

    gK = R^2 - x^2 - y^2;
    gX = 1 - x^2 - y^2;
    gS = s*(1-s);

    [V,cV] = polynomial([x y s], degV);

    degMult = min(max(degV-2,0),2);

    [m1K,c1K] = polynomial([x y s], degMult);
    [m1S,c1S] = polynomial([x y s], degMult);

    [m2K,c2K] = polynomial([x y s], degMult);
    [m2S,c2S] = polynomial([x y s], degMult);

    [m3X,c3X] = polynomial([x y], degMult);

    Vdot = jacobian(V,x)*(T*f1) ...
         + jacobian(V,y)*(T*f2) ...
         + jacobian(V,s);

    D1 = -Vdot;
    D2 = V - Phi;
    D3 = C - replace(V,s,0);

    constraints = [
        sos(m1K), ...
        sos(m1S), ...
        sos(m2K), ...
        sos(m2S), ...
        sos(m3X), ...
        sos(D1 - m1K*gK - m1S*gS), ...
        sos(D2 - m2K*gK - m2S*gS), ...
        sos(D3 - m3X*gX)
        ];

    params = [cV; c1K; c1S; c2K; c2S; c3X; C];

    opts = sdpsettings('solver','mosek','verbose',verbose);

    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-7;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-7;
    opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-7;

    sol = solvesos(constraints, C, opts, params);

    Cval = value(C);

    if nargout > 2

        assign(params, value(params));

        aux.vars = [x y s];
        aux.Vnum = clean(replace(V, params, value(params)), 1e-9);
        aux.D1num = clean(replace(D1, params, value(params)), 1e-9);
        aux.D2num = clean(replace(D2, params, value(params)), 1e-9);
        aux.C = Cval;
    end
end

function val = exact_extreme(a,T)

    ts = candidate_times(a,T);
    vals = exp(-2*ts).*(1 + a^2*ts.^2);

    val = max(vals);
end

function tstar = maximizing_time(a,T)

    ts = candidate_times(a,T);
    vals = exp(-2*ts).*(1 + a^2*ts.^2);

    [~,idx] = max(vals);

    tstar = ts(idx);
end

function ts = candidate_times(a,T)

    ts = [0,T];

    if a > 2

        disc = 1 - 4/a^2;

        tcrit1 = (1 - sqrt(disc))/2;
        tcrit2 = (1 + sqrt(disc))/2;

        ts = [ts,tcrit1,tcrit2];
    end

    ts = unique(ts(ts >= 0 & ts <= T));
end

function val = eval_poly(p, vars, xyz)

    val = double(replace(p, vars, xyz));
end

function plot_defect_field(X,Y,Dvals,theta,xstar,ystar,idxStar, ...
    x0s,y0s,S,fileName,saveFigures)

    fig = figure;
    set(fig,'Color','w','Units','centimeters','Position',[2 2 11 10]);

    hold on
    box on

    contourf(X,Y,log10(Dvals),35,'LineColor','none');
    colorbar
    colormap(parula)

    plot(cos(theta), sin(theta), '--', ...
        'Color', S.black, ...
        'LineWidth', 2);

    plot(xstar,ystar, ...
        'Color',S.red, ...
        'LineWidth',5.0);

    arrow_ids = round([0.35 0.62 0.88]*length(xstar));

    for jj = arrow_ids
        add_arrow_head_on_curve(xstar,ystar,jj,0.07,S.red);
    end

    scatter(xstar(idxStar), ystar(idxStar), 120, S.red, 'filled');

    scatter(x0s,y0s,120,S.red,'filled', ...
        'MarkerEdgeColor',0.55*S.red, ...
        'LineWidth',0.8);

    scatter(0,0,100,S.black,'filled');

    set(gca,'FontSize',42);

    xlabel('$x$','Interpreter','latex','FontSize',42,'FontWeight','bold');
    ylabel('$y$','Interpreter','latex','FontSize',42,'FontWeight','bold');

    if saveFigures
        save_pdf_figure(fig,fileName);
    end
end

function add_arrow_head_on_curve(xtraj,ytraj,idx,sz,col)

    idx = max(3,min(idx,length(xtraj)-3));

    dx = xtraj(idx+2)-xtraj(idx-2);
    dy = ytraj(idx+2)-ytraj(idx-2);

    ang = atan2(dy,dx);

    tri = sz*[1 0; -0.65 0.38; -0.65 -0.38];

    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    triR = (R*tri')';

    patch(xtraj(idx)+triR(:,1), ytraj(idx)+triR(:,2), col, ...
        'EdgeColor','none', ...
        'FaceAlpha',1);
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