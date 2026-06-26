%% rossler_poincare_periodic_orbits.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Data-Driven Auxiliary Function Methods
% Example: Data-driven moment constraints for Rössler periodic orbits
%
% This script identifies unstable periodic orbits of the Rössler system
% using data from a Poincaré return map.
%
% The script:
%
%   1. simulates a long trajectory of the Rössler system,
%   2. extracts crossings of the Poincaré section x_1 = 0 with x_1' > 0,
%   3. learns a Chebyshev EDMD approximation of the return map,
%   4. imposes data-driven invariant moment constraints,
%   5. solves random moment-SDPs to extract atomic invariant measures, and
%   6. plots the recovered periodic points on iterates of the map.
%
% Pressing Run should reproduce the numerical output and save:
%
%   rossler_poincare_data.pdf
%   rossler_poincare_scaled_data.pdf
%   rossler_poincare_period2.pdf
%   rossler_poincare_period3.pdf
%   rossler_poincare_period4.pdf
%   rossler_poincare_periodic_orbits.mat
%
% Requirements:
%   - MATLAB
%   - YALMIP
%   - MOSEK
%   - Chebfun
%
% -------------------------------------------------------------------------

clear; clc; close all;
yalmip('clear');
format long;
rng(1);

%% User options

a = 0.1;
b = 0.1;
c = 18;

numSteps = 1e6;
dt = 0.005;
x0 = [0; -15; 0];

chebDegreePsi = 20;
chebDegreeTheta = 4*chebDegreePsi;

edmdTol = 0;

numRandomObjectives = 1000;
maxStoredPeriod = 100;

saveFigures = true;
saveData = true;

%% Book plotting style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.green  = [0, 140, 90]/255;
S.orange = [230, 120, 20]/255;
S.purple = [170, 90, 160]/255;
S.black  = [0, 0, 0];
S.gray   = [0.50, 0.50, 0.50];

%% Simulate Rössler trajectory

fprintf('\n============================================================\n');
fprintf('Rössler Poincaré periodic orbits from moment constraints\n');
fprintf('============================================================\n\n');

fprintf('Parameters:\n');
fprintf('  a = %.6g\n', a);
fprintf('  b = %.6g\n', b);
fprintf('  c = %.6g\n', c);
fprintf('numSteps = %d, dt = %.6g\n\n', numSteps, dt);

rossler = @(~,x) [
    -x(2) - x(3);
     x(1) + a*x(2);
     b + x(3)*(x(1) - c)
];

odeopts = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
tspan = (0:numSteps)*dt;

fprintf('Integrating Rössler trajectory...\n');

[~,xdat] = ode45(rossler,tspan,x0,odeopts);

%% Collect Poincaré section data

fprintf('Collecting Poincaré section crossings...\n');

sectionY = [];
sectionZ = [];
sectionT = [];

for ind = 1:numSteps

    if xdat(ind,1) < 0 && xdat(ind+1,1) >= 0

        alpha = -xdat(ind,1)/(xdat(ind+1,1) - xdat(ind,1));

        xcross = xdat(ind,:) ...
            + alpha*(xdat(ind+1,:) - xdat(ind,:));

        tcross = tspan(ind) ...
            + alpha*(tspan(ind+1) - tspan(ind));

        fcross = rossler(tcross,xcross(:));

        if fcross(1) > 0
            sectionY(end+1,1) = xcross(2); %#ok<SAGROW>
            sectionZ(end+1,1) = xcross(3); %#ok<SAGROW>
            sectionT(end+1,1) = tcross;    %#ok<SAGROW>
        end
    end
end

x = sectionY;
numSectionPoints = numel(x);

fprintf('Collected %d Poincaré section crossings.\n',numSectionPoints);
fprintf('Max |z| on section: %.3e\n',max(abs(sectionZ)));

if numSectionPoints < 2
    error('Not enough Poincaré section crossings were collected.');
end

%% Scale section coordinate to [-1,1]

xmin = min(x);
xmax = max(x);

xs = 2*(x - xmin)/(xmax - xmin) - 1;

fprintf('\nSection coordinate scaling:\n');
fprintf('  xmin = %.12g\n', xmin);
fprintf('  xmax = %.12g\n', xmax);

%% Figure 1: return map data

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 11 9]);

hold on
box on

plot(x,x,'Color',S.gray,'LineWidth',2);
plot(x(1:end-1),x(2:end),'k.','MarkerSize',18);

xlabel('$x_{2,n}$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$x_{2,n+1}$','Interpreter','latex','FontSize',24,'FontWeight','bold');

axis tight
set(gca,'FontSize',24,'TickLabelInterpreter','latex');
grid on

if saveFigures
    save_pdf_figure(fig1,'rossler_poincare_data.pdf');
end

%% Figure 2: scaled return map data

fig2 = figure;
set(fig2,'Color','w','Units','centimeters','Position',[2 2 11 9]);

hold on
box on

plot(xs(1:end-1),xs(2:end),'k.','MarkerSize',18);

xlabel('$\tilde{x}_{2,n}$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$\tilde{x}_{2,n+1}$','Interpreter','latex','FontSize',24,'FontWeight','bold');

axis tight
set(gca,'FontSize',24,'TickLabelInterpreter','latex');
grid on

if saveFigures
    save_pdf_figure(fig2,'rossler_poincare_scaled_data.pdf');
end

%% Build Chebyshev EDMD approximation and discrete-time generator

fprintf('\nBuilding Chebyshev EDMD approximation...\n');

if mod(chebDegreeTheta,2) ~= 0
    error('chebDegreeTheta must be even.');
end

Tpsi   = chebpoly(0:chebDegreePsi,[-1,1]);
Ttheta = chebpoly(0:chebDegreeTheta,[-1,1]);

Phi = Ttheta(xs(1:end-1))';
Psi = Tpsi(xs(2:end))';

K = edmd_with_thresholding(Phi,Psi,edmdTol);

% Discrete-time Poincaré map with tau = 1:
%     L = K - I,
% where I is the rectangular identity embedding.
L = K - eye(size(K));

%% Moment and localizing matrices in Chebyshev basis

fprintf('Constructing Chebyshev moment matrices...\n');

y = sdpvar(chebDegreeTheta,1);

A = chebsdp_1d(chebDegreeTheta/2);
M0 = reshape(A*[1;y], ...
    [chebDegreeTheta/2+1, chebDegreeTheta/2+1]);

Bmat = chebsdp_1d_locball(chebDegreeTheta/2);
M1 = reshape(Bmat*[1;y], ...
    [chebDegreeTheta/2, chebDegreeTheta/2]);

%% Solve random moment SDPs

fprintf('Solving %d random moment SDPs...\n',numRandomObjectives);

upos = zeros(numRandomObjectives,maxStoredPeriod+1);

solverOpts = sdpsettings('solver','mosek','verbose',0);

constraints = [
    L*[1;y] == 0, ...
    M0 >= 0, ...
    M1 >= 0
];

for ind = 1:numRandomObjectives

    objective = dot(10*rand(1,length(y(1:5))) - 5, y(1:5));

    sol = optimize(constraints,objective,solverOpts);

    if sol.problem ~= 0
        fprintf('Iteration %d: solver status %d (%s)\n', ...
            ind,sol.problem,sol.info);
        continue
    end

    M0val = value(M0);

    if any(isnan(M0val(:)))
        fprintf('Iteration %d: moment matrix contains NaNs.\n',ind);
        continue
    end

    atomsScaled = cheb_extract_minimizers(M0val,0:chebDegreeTheta/2);
    atomsScaled = sort(real(atomsScaled(:)))';

    atomsScaled = atomsScaled(atomsScaled >= -1-1e-8 ...
                            & atomsScaled <=  1+1e-8);

    atoms = unscale_section_coordinate(atomsScaled,xmin,xmax);

    period = numel(atoms);
    keep = min(period,maxStoredPeriod);

    if keep > 0
        upos(ind,1:keep+1) = [period, atoms(1:keep)];
    end

    fprintf('Iteration %4d: extracted %d section points.\n',ind,period);
end

%% Periodic orbit points used in the book figures

period1Points = -22.9179234809962;

period2Points = [
    -18.0313858364717
    -24.4815932442540
];

period3A = [
    -25.7885341999335
    -14.0390147030195
    -19.1857534188735
];

period3B = [
    -16.1726475968098
    -25.0761646600292
    -22.0681328549753
];

period4Points = [
    -16.9994950231043
    -22.2068733571595
    -23.1419355614561
    -24.8175866399386
];

%% Plot selected iterates

plot_iterate_figure( ...
    x, 2, ...
    {period1Points, period2Points}, ...
    {S.red, S.blue}, ...
    'rossler_poincare_period2.pdf', ...
    saveFigures);

plot_iterate_figure( ...
    x, 3, ...
    {period1Points, period3A, period3B}, ...
    {S.red, S.orange, S.purple}, ...
    'rossler_poincare_period3.pdf', ...
    saveFigures);

plot_iterate_figure( ...
    x, 4, ...
    {period1Points, period2Points, period4Points}, ...
    {S.red, S.blue, S.green}, ...
    'rossler_poincare_period4.pdf', ...
    saveFigures);

%% Save output

if saveData
    save('rossler_poincare_periodic_orbits.mat', ...
        'upos','x','xs','sectionZ','sectionT', ...
        'xmin','xmax','K','L', ...
        'chebDegreePsi','chebDegreeTheta','edmdTol', ...
        'numRandomObjectives','maxStoredPeriod', ...
        'period1Points','period2Points','period3A','period3B','period4Points');
    fprintf('\nSaved data: rossler_poincare_periodic_orbits.mat\n');
end

fprintf('\nFinished Rössler Poincaré periodic-orbit example.\n');

%% Local helper functions

function val = unscale_section_coordinate(xs,xmin,xmax)

    val = 0.5*(xs + 1)*(xmax - xmin) + xmin;
end

function plot_iterate_figure(x,iterate,pointSets,colors,fileName,saveFigure)

    fig = figure;
    set(fig,'Color','w','Units','centimeters','Position',[2 2 11 9]);

    hold on
    box on

    plot(x,x,'Color',[0.5 0.5 0.5],'LineWidth',2);
    plot(x(1:end-iterate),x(1+iterate:end),'k.','MarkerSize',18);

    for j = 1:numel(pointSets)

        pts = pointSets{j};

        plot(pts,pts,'.', ...
            'Color',colors{j}, ...
            'MarkerSize',45);
    end

    xlabel(sprintf('$x_{2,n}$'), ...
        'Interpreter','latex','FontSize',24,'FontWeight','bold');

    ylabel(sprintf('$x_{2,n+%d}$',iterate), ...
        'Interpreter','latex','FontSize',24,'FontWeight','bold');

    axis tight
    set(gca,'FontSize',24,'TickLabelInterpreter','latex');
    grid on

    if saveFigure
        save_pdf_figure(fig,fileName);
    end
end

function K = edmd_with_thresholding(Phi,Psi,tol)

    if nargin < 3
        tol = 0;
    end

    [U,S,V] = svd(Phi,'econ');
    s = diag(S);

    if isempty(s)
        error('Phi has no singular values.');
    end

    if tol > 0
        keep = s > tol*max(s);
    else
        keep = s > eps(max(size(Phi)))*max(s);
    end

    if ~any(keep)
        error('All singular values were discarded. Decrease tol.');
    end

    PhiPinv = V(:,keep)*diag(1./s(keep))*U(:,keep)';
    K = Psi*PhiPinv;
end

function x = cheb_extract_minimizers(momentMatrices,gramMonomials)

    tol = 1e-12;
    dropTol = 1e-1;

    if iscell(momentMatrices)

        x = cell(size(momentMatrices));

        for kk = 1:length(momentMatrices)

            Mk = momentMatrices{kk};
            cleanTol = tol*max(abs(Mk(:)));

            Mk = clean(Mk,cleanTol);

            [U,S] = svd(Mk);

            [Sdiag,pos] = sort(diag(S),'descend');
            U = U(:,pos);

            drop = Sdiag(2:end)./(eps + Sdiag(1:end-1));
            drop = find(drop < dropTol,1,'first');

            if ~isempty(drop)
                rankM = drop;
            else
                rankM = length(Sdiag);
            end

            if rankM == 0
                x{kk} = [];
                continue
            end

            U = U(:,1:rankM)*diag(sqrt(Sdiag(1:rankM)));

            [U,basis] = column_echelon_form(U,1e-6);

            beta = gramMonomials{kk}(basis);
            nmons = length(gramMonomials{kk});

            I = [];
            J = [];

            for jj = 1:rankM

                powers = [beta(jj)+1; abs(beta(jj)-1)];

                [~,pos2] = ismember(powers(:),gramMonomials{kk}(:));

                I = [I; jj; jj]; %#ok<AGROW>
                J = [J; pos2(:)]; %#ok<AGROW>
            end

            multiplicationMatrix = sparse( ...
                I,J,0.5*ones(size(I)),rankM,nmons);

            Nmat = multiplicationMatrix*U;

            [Q,~] = schur(full(Nmat));

            if isempty(Q)
                x{kk} = [];
            else
                x{kk} = zeros(1,rankM);

                for ii = 1:rankM
                    x{kk}(ii) = Q(:,ii)'*Nmat*Q(:,ii);
                end
            end
        end

    else

        xcell = cheb_extract_minimizers({momentMatrices},{gramMonomials});
        x = xcell{1};
    end
end

function At = chebsdp_1d(d)

    [i,j] = meshgrid(1:d+1,1:d+1);

    s = sub2ind([d+1,d+1],i,j);

    i = i - 1;
    j = j - 1;

    ipj = i + j;
    imj = abs(i - j);

    row = [];
    col = [];
    val = [];

    for kk = 0:2*d

        idx = [s(ipj == kk); s(imj == kk)];

        row = [row; idx]; %#ok<AGROW>
        col = [col; kk*ones(size(idx)) + 1]; %#ok<AGROW>
        val = [val; 0.5*ones(size(idx))]; %#ok<AGROW>
    end

    At = sparse(row,col,val,(d+1)^2,2*d+1);
end

function Bt = chebsdp_1d_locball(d)

    At = chebsdp_1d(d-1);
    Bt = At;
    Bt(:,end+1:end+2) = 0;

    T2 = chebpoly(2);

    for ii = 1:size(Bt,2)-2

        Ti = chebpoly(ii-1);
        TiT2 = chebcoeffs(T2.*Ti);
        TiT2 = TiT2(:).';

        TiT2(abs(TiT2) < 1e-14) = 0;

        tmp = At(:,ii).*TiT2;

        Bt(:,1:ii+2) = Bt(:,1:ii+2) - tmp;
    end
end

function [A,basis] = column_echelon_form(A,tol)

    [n,m] = size(A);

    i = 1;
    j = 1;
    basis = [];

    while (i <= m) && (j <= n)

        [pivot,k] = max(abs(A(j,i:m)));
        k = k + i - 1;

        if pivot <= tol

            A(j,i:m) = zeros(1,m-i+1);
            j = j + 1;

        else

            basis = [basis j]; %#ok<AGROW>
            A(j:n,[i k]) = A(j:n,[k i]);

            found = false;

            while ~found

                if abs(A(j,i)) < tol*max(abs(A(:,i)))
                    j = j + 1;
                    found = (j == n);
                else
                    found = true;
                end
            end

            if j <= n

                A(j:n,i) = A(j:n,i)/A(j,i);

                for kk = [1:i-1 i+1:m]
                    A(j:n,kk) = A(j:n,kk) - A(j,kk)*A(j:n,i);
                end

                i = i + 1;
                j = j + 1;
            end
        end
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