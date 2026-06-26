%% cylinder_generate_pod_data.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Data-Driven Auxiliary Function Methods
% Example: POD data generation for cylinder wake ROM learning
%
% Simulates 2D incompressible flow past a circular cylinder at Re = 100,
% stores vorticity snapshots, computes a POD basis from the settled wake,
% and projects all snapshots onto the first POD modes.
%
% Output:
%   cylinder_Re100_POD10_projected_all.mat
%   cylinder_wake_snapshots.pdf
%   cylinder_wake_pod_modes.pdf
%   cylinder_Re100_vorticity_movie.mp4
%
% Requirements:
%   - MATLAB
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% User options

Re = 100;

Lx = 18;
Ly = 5;
D  = 1;

Umax  = 1.0;
Umean = (2/3)*Umax;
nu    = Umean*D/Re;

xc = 4.0;
yc = Ly/2;
Rcyl = D/2;

Nx = 720;
Ny = 200;

dt = 0.001;
Tfinal = 200;
snapEvery = 100;

rPOD = 10;
tPODStart = 30;

plotTimes = [20 40];
nModesPlot = 10;

saveFigures = true;
saveData = true;
saveMovie = true;
playMovieAtEnd = false;

%% Snapshot grid

Nx_snap = 4957;
Ny_snap = 32;

x_snap = linspace(0,Lx,Nx_snap);
y_snap = linspace(0,Ly,Ny_snap);

[Xs,Ys] = meshgrid(x_snap,y_snap);

%% Smooth plotting grid

Nx_plot = 1200;
Ny_plot = 350;

x_plot = linspace(0,Lx,Nx_plot);
y_plot = linspace(0,Ly,Ny_plot);

[Xp,Yp] = meshgrid(x_plot,y_plot);

%% Computational grid

x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);

dx = x(2)-x(1);
dy = y(2)-y(1);

[X,Y] = meshgrid(x,y);

solid = (X-xc).^2 + (Y-yc).^2 <= Rcyl^2;

%% Time stepping

nSteps = round(Tfinal/dt);

uin = Umax*4*y.*(Ly-y)/Ly^2;
uin = uin(:);

uvel = repmat(uin,1,Nx);
vvel = zeros(Ny,Nx);

rng(1);
vvel = vvel + 1e-2*randn(size(vvel));

vvel(:,1) = 0;
uvel(solid) = 0;
vvel(solid) = 0;

p = zeros(Ny,Nx);

snapshots = [];
snapTimes = [];
smoothSnapshots = [];

movieFrames = struct('cdata',{},'colormap',{});
frameCount = 0;

figMovie = figure('Color','w');
set(figMovie,'Renderer','opengl');

fprintf('\n============================================================\n');
fprintf('Cylinder wake POD data generation\n');
fprintf('============================================================\n\n');

fprintf('Re = %g\n',Re);
fprintf('Grid: Nx = %d, Ny = %d\n',Nx,Ny);
fprintf('dt = %.4g, Tfinal = %.4g\n',dt,Tfinal);
fprintf('Snapshot interval = %.4g\n\n',snapEvery*dt);

%% Main simulation loop

for n = 1:nSteps

    t = n*dt;

    [ux,uy] = grad_field(uvel,dx,dy);
    [vx,vy] = grad_field(vvel,dx,dy);

    Lu = lap_field(uvel,dx,dy);
    Lv = lap_field(vvel,dx,dy);

    ustar = uvel + dt*(-uvel.*ux - vvel.*uy + nu*Lu);
    vstar = vvel + dt*(-uvel.*vx - vvel.*vy + nu*Lv);

    [ustar,vstar] = apply_velocity_bc(ustar,vstar,uin,solid);

    divstar = div_field(ustar,vstar,dx,dy);
    rhs = divstar/dt;

    p = poisson_pressure(p,rhs,dx,dy,solid,200);

    [px,py] = grad_field(p,dx,dy);

    uvel = ustar - dt*px;
    vvel = vstar - dt*py;

    [uvel,vvel] = apply_velocity_bc(uvel,vvel,uin,solid);

    if t < 15
        wake = exp(-((X-(xc+1.0)).^2/0.5^2 + (Y-yc).^2/0.35^2));
        vvel = vvel + dt*0.35*wake.*sin(2*pi*0.25*t);

        vvel(:,1) = 0;
        uvel(solid) = 0;
        vvel(solid) = 0;
    end

    if mod(n,snapEvery) == 0

        omega = vorticity_field(uvel,vvel,dx,dy);
        omega(solid) = NaN;

        omegaSnap = interp2(X,Y,omega,Xs,Ys,'linear',0);

        snapshots(:,end+1) = omegaSnap(:); %#ok<SAGROW>
        snapTimes(end+1) = t;              %#ok<SAGROW>

        omegaPlotField = omega;
        omegaPlotField(~isfinite(omegaPlotField)) = 0;
        omegaPlotField(solid) = 0;

        omegaPlot = interp2(X,Y,omegaPlotField,Xp,Yp,'linear',0);
        smoothSnapshots(:,end+1) = omegaPlot(:); %#ok<SAGROW>

        clf(figMovie);

        imagesc(x_plot,y_plot,omegaPlot);
        set(gca,'YDir','normal');
        axis equal tight;

        colormap(parula);
        colorbar;
        clim([-1 1]);

        hold on

        th = linspace(0,2*pi,400);
        xcyl = xc + Rcyl*cos(th);
        ycyl = yc + Rcyl*sin(th);

        fill(xcyl,ycyl,[0.95 0.95 0.95], ...
            'EdgeColor','k', ...
            'LineWidth',1.5);

        hold off

        xlabel('x');
        ylabel('y');
        title(sprintf('Cylinder wake, Re = %d, t = %.2f',Re,t));

        drawnow;

        frameCount = frameCount + 1;
        movieFrames(frameCount) = getframe(figMovie); %#ok<SAGROW>

        fprintf('t = %.3f, snapshot size = %d values\n',t,numel(omegaSnap));
    end
end

%% Plot selected vorticity snapshots

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 16 5]);

tiledlayout(1,numel(plotTimes), ...
    'TileSpacing','compact', ...
    'Padding','loose');

for jj = 1:numel(plotTimes)

    [~,idx] = min(abs(snapTimes - plotTimes(jj)));
    tplot = snapTimes(idx);

    omegaPlot = reshape(smoothSnapshots(:,idx),Ny_plot,Nx_plot);

    nexttile;

    imagesc(x_plot,y_plot,omegaPlot);
    set(gca,'YDir','normal');

    axis equal tight;

    colormap(parula);
    clim([-1 1]);

    hold on

    th = linspace(0,2*pi,400);
    xcyl = xc + (Rcyl+0.15)*cos(th);
    ycyl = yc + (Rcyl+0.15)*sin(th);

    fill(xcyl,ycyl,[0.95 0.95 0.95], ...
        'EdgeColor','k', ...
        'LineWidth',1.2);

    hold off
    box on

    xlabel('$x$','Interpreter','latex','FontSize',24);
    ylabel('$y$','Interpreter','latex','FontSize',24);

    title(sprintf('$t = %.0f$',tplot), ...
        'Interpreter','latex');

    set(gca,'FontSize',18,'TickLabelInterpreter','latex');
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontSize = 16;
cb.Position(1) = cb.Position(1) - 0.01;

if saveFigures
    export_pdf(fig1,'cylinder_wake_snapshots.pdf');
end

%% POD basis from post-transient data

Xall = snapshots;
Xall(~isfinite(Xall)) = 0;

idxPOD = snapTimes >= tPODStart;
XpodBasis = Xall(:,idxPOD);

omegaMean = mean(XpodBasis,2);

XflucBasis = XpodBasis - omegaMean;

[Phi,Svals,~] = svd(XflucBasis,'econ');

Phi_r = Phi(:,1:rPOD);

XallFluc = Xall - omegaMean;
A = (Phi_r.'*XallFluc).';

singVals = diag(Svals);
energyFrac = singVals.^2/sum(singVals.^2);
cumEnergy = cumsum(energyFrac);

fprintf('\nPOD basis computed using snapshots with t >= %.2f\n',tPODStart);
fprintf('Projected all %d snapshots onto first %d POD modes\n',numel(snapTimes),rPOD);
fprintf('Settled-wake POD energy captured by first %d modes: %.6f\n', ...
    rPOD,cumEnergy(rPOD));

dataFileOut = sprintf('cylinder_Re%d_POD10_projected_all.mat',Re);

if saveData
    save(dataFileOut, ...
        'Phi_r','A','omegaMean','singVals','energyFrac','cumEnergy', ...
        'snapTimes','tPODStart','x_snap','y_snap','rPOD','Re', ...
        'Lx','Ly','D','Umax','Umean','nu','xc','yc', ...
        '-v7.3');

    fprintf('Saved data: %s\n',dataFileOut);
end

%% Plot POD modes

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 16 15]);

tiledlayout(5,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for j = 1:nModesPlot

    modeSnap = reshape(Phi_r(:,j),Ny_snap,Nx_snap);
    modePlot = interp2(Xs,Ys,modeSnap,Xp,Yp,'linear',0);

    nexttile;

    imagesc(x_plot,y_plot,modePlot);
    set(gca,'YDir','normal');

    axis equal tight;
    colormap(parula);

    cmax = max(abs(modePlot(:)));
    clim([-cmax cmax]);

    hold on

    th = linspace(0,2*pi,400);
    xcyl = xc + (Rcyl+0.15)*cos(th);
    ycyl = yc + (Rcyl+0.15)*sin(th);

    fill(xcyl,ycyl,[0.95 0.95 0.95], ...
        'EdgeColor','k', ...
        'LineWidth',1.2);

    hold off

    xlabel('$x$','Interpreter','latex','FontSize',24);
    ylabel('$y$','Interpreter','latex','FontSize',24);

    colorbar
    box on

    title(sprintf('POD mode %d',j), ...
        'Interpreter','latex');

    set(gca,'FontSize',18,'TickLabelInterpreter','latex');
end

if saveFigures
    export_pdf(fig2,'cylinder_wake_pod_modes.pdf');
end

%% Movie

if playMovieAtEnd
    figure('Color','w');
    movie(gcf,movieFrames,1,20);
end

if saveMovie
    movieFileOut = sprintf('cylinder_Re%d_vorticity_movie.mp4',Re);

    vobj = VideoWriter(movieFileOut,'MPEG-4');
    open(vobj);

    for k = 1:numel(movieFrames)
        writeVideo(vobj,movieFrames(k));
    end

    close(vobj);

    fprintf('Saved movie: %s\n',movieFileOut);
end

fprintf('\nFinished cylinder POD-data generation.\n');

%% Local functions

function [fx,fy] = grad_field(f,dx,dy)

    fx = zeros(size(f));
    fy = zeros(size(f));

    fx(:,2:end-1) = (f(:,3:end)-f(:,1:end-2))/(2*dx);
    fx(:,1)       = (f(:,2)-f(:,1))/dx;
    fx(:,end)     = (f(:,end)-f(:,end-1))/dx;

    fy(2:end-1,:) = (f(3:end,:)-f(1:end-2,:))/(2*dy);
    fy(1,:)       = (f(2,:)-f(1,:))/dy;
    fy(end,:)     = (f(end,:)-f(end-1,:))/dy;
end

function L = lap_field(f,dx,dy)

    L = zeros(size(f));

    L(2:end-1,2:end-1) = ...
        (f(2:end-1,3:end)-2*f(2:end-1,2:end-1)+f(2:end-1,1:end-2))/dx^2 + ...
        (f(3:end,2:end-1)-2*f(2:end-1,2:end-1)+f(1:end-2,2:end-1))/dy^2;
end

function d = div_field(u,v,dx,dy)

    [ux,~] = grad_field(u,dx,dy);
    [~,vy] = grad_field(v,dx,dy);

    d = ux + vy;
end

function w = vorticity_field(u,v,dx,dy)

    [~,uy] = grad_field(u,dx,dy);
    [vx,~] = grad_field(v,dx,dy);

    w = vx - uy;
end

function [u,v] = apply_velocity_bc(u,v,uin,solid)

    [Ny,Nx] = size(u);

    u(:,1) = uin;
    v(:,1) = 0;

    u(:,Nx) = u(:,Nx-1);
    v(:,Nx) = v(:,Nx-1);

    u(1,:)  = 0;
    v(1,:)  = 0;

    u(Ny,:) = 0;
    v(Ny,:) = 0;

    u(solid) = 0;
    v(solid) = 0;
end

function p = poisson_pressure(p,rhs,dx,dy,solid,nIter)

    beta = dx^2*dy^2/(2*(dx^2+dy^2));

    for k = 1:nIter %#ok<NASGU>

        pold = p;

        p(2:end-1,2:end-1) = beta*( ...
            (pold(2:end-1,3:end)+pold(2:end-1,1:end-2))/dx^2 + ...
            (pold(3:end,2:end-1)+pold(1:end-2,2:end-1))/dy^2 - ...
            rhs(2:end-1,2:end-1));

        p(:,1)   = p(:,2);
        p(:,end) = 0;

        p(1,:)   = p(2,:);
        p(end,:) = p(end-1,:);

        p(solid) = 0;
    end
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