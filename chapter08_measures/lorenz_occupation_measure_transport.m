%% lorenz_occupation_measure_transport.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Invariant Measures, Ergodic Optimization, and Duality
% Example: Occupation-measure transport for the Lorenz system
%
% This script computes finite-time terminal measures for the Lorenz system
% using a time-state occupation-measure moment SDP.
%
% The Lorenz system is first rescaled to a unit cube containing a numerically
% observed attractor. The initial measure is a uniform probability measure
% on a small ball in one lobe. For each terminal time T, the script solves
% the Liouville moment constraints
%
% Time is normalized by s=t/T, so the implemented constraints are
%
%     int (partial_s v + T*grad_x(v).f) deta_bar
%         = int v(1,x) dmu_T - int v(0,x) dmu_0.
%
% and reconstructs a polynomial density for the terminal measure.
%
% Pressing Run should reproduce the numerical output and save:
%
%     lorenz_occ_density_T_0p5.pdf
%     lorenz_occ_density_T_1.pdf
%     lorenz_occ_density_T_5.pdf
%     lorenz_occ_density_T_10.pdf
%
% Requirements:
%   - MATLAB
%   - YALMIP
%   - MOSEK
%
% -------------------------------------------------------------------------

clear; clc; close all;
yalmip('clear');
rng(1);

%% User options

sigma = 10;
rho   = 28;
beta  = 8/3;

momentDegree = 10;
terminalTimes = [0.5 1 5 10];

Tattr = 160;
dtAttr = 0.001;

initialCenterScaled = [0.746040; 0.331799; 0.668536];
initialRadiusScaled = 0.025;

numInitialPlotPoints = 1200;
numValidationPoints = 250;

solver_name = 'mosek';
verbose_solver = 1;

saveFigures = true;

%% Solver options

optsSDP = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver);

%% Build a scaled Lorenz attractor for visualization and support

fprintf('\n============================================================\n');
fprintf('Lorenz occupation-measure transport\n');
fprintf('============================================================\n\n');

fprintf('Lorenz parameters:\n');
fprintf('  sigma = %.6g\n', sigma);
fprintf('  rho   = %.6g\n', rho);
fprintf('  beta  = %.6g\n', beta);
fprintf('Moment degree d = %d\n\n', momentDegree);

fprintf('Simulating Lorenz attractor for scaling and plotting...\n');

optsODE = odeset('RelTol',1e-9,'AbsTol',1e-11);
tspan = 0:dtAttr:Tattr;

[~, Xraw] = ode45(@(t,x) lorenz_rhs(t,x,sigma,rho,beta), ...
    tspan, [1;1;1], optsODE);

Xraw = Xraw(round(end/3):end,:);

xmin = min(Xraw,[],1);
xmax = max(Xraw,[],1);

pad = 0.12*(xmax - xmin);

xmin = xmin - pad;
xmax = xmax + pad;

center = (xmax + xmin)/2;
scale  = (xmax - xmin)/2;

Yattr = (Xraw - center)./scale;

fprintf('\nScaling to unit cube K = [-1,1]^3:\n');
fprintf('  center = [% .8f % .8f % .8f]\n', center);
fprintf('  scale  = [% .8f % .8f % .8f]\n', scale);
fprintf('  max |Yattr_i| = [% .4f % .4f % .4f]\n\n', ...
    max(abs(Yattr),[],1));

fprintf('Initial ball in scaled variables:\n');
fprintf('  center = [% .6f % .6f % .6f]\n', initialCenterScaled);
fprintf('  radius = %.6f\n\n', initialRadiusScaled);

%% Polynomial data

vectorFieldTerms = scaled_lorenz_polynomial_terms( ...
    center(:), scale(:), sigma, rho, beta);

basisState = monomial_basis_3d(momentDegree);
numStateMoments = size(basisState,1);
idxState = make_index_map_3d(basisState);

basisOcc = monomial_basis_4d(momentDegree);
numOccMoments = size(basisOcc,1);
idxOcc = make_index_map_4d(basisOcc);

%% Prescribed initial moments

initialMoments = zeros(numStateMoments,1);

for k = 1:numStateMoments
    initialMoments(k) = shifted_ball_moment( ...
        basisState(k,:), initialCenterScaled(:).', initialRadiusScaled);
end

Y0validation = sample_ball(numValidationPoints, ...
    initialCenterScaled(:).', initialRadiusScaled);
X0validation = Y0validation.*scale + center;

%% Solve occupation-measure SDP for each terminal time

terminalMoments = cell(numel(terminalTimes),1);

for kk = 1:numel(terminalTimes)

    T = terminalTimes(kk);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Solving occupation-measure SDP for T = %g\n', T);
    fprintf('------------------------------------------------------------\n');

    yOcc  = sdpvar(numOccMoments,1);
    yTerm = sdpvar(numStateMoments,1);

    constraints = [];

    %% Liouville constraints

    constraints = [constraints, build_liouville_constraints( ...
        yOcc, yTerm, initialMoments, basisOcc, idxOcc, idxState, ...
        vectorFieldTerms, momentDegree, T)];

    %% Mass constraints

    constraints = [
        constraints, ...
        moment_expr_4d(yOcc,idxOcc,0,0,0,0) == 1, ...
        moment_expr_3d(yTerm,idxState,0,0,0) == 1
        ];

    %% Support constraints on K = [-1,1]^3

    basisMomentOcc = monomial_basis_4d(floor(momentDegree/2));
    basisMomentState = monomial_basis_3d(floor(momentDegree/2));

    constraints = [
        constraints, ...
        build_moment_matrix_4d(yOcc,idxOcc,basisMomentOcc) >= 0, ...
        build_moment_matrix_3d(yTerm,idxState,basisMomentState) >= 0
        ];

    basisLocalizingOcc = monomial_basis_4d(floor((momentDegree-2)/2));
    basisLocalizingState = monomial_basis_3d(floor((momentDegree-2)/2));

    constraints = [constraints, ...
        build_time_localizing_matrix_4d( ...
            yOcc,idxOcc,basisLocalizingOcc) >= 0];

    for coord = 1:3
        constraints = [
            constraints, ...
            build_cube_localizing_matrix_4d( ...
                yOcc,idxOcc,basisLocalizingOcc,coord) >= 0, ...
            build_cube_localizing_matrix_3d( ...
                yTerm,idxState,basisLocalizingState,coord) >= 0
            ]; %#ok<AGROW>
    end

    %% Small regularization for reproducible feasible terminal moments

    objective = 1e-9*(yOcc.'*yOcc + yTerm.'*yTerm);

    diagnostics = optimize(constraints,objective,optsSDP);

    fprintf('  status: %s\n', diagnostics.info);

    if diagnostics.problem ~= 0
        warning('SDP at T = %g did not solve cleanly.', T);
    end

    yTermValue = value(yTerm);
    terminalMoments{kk} = yTermValue;

    terminalMass = yTermValue(get_idx_3d(idxState,0,0,0));

    terminalMean = [
        yTermValue(get_idx_3d(idxState,1,0,0));
        yTermValue(get_idx_3d(idxState,0,1,0));
        yTermValue(get_idx_3d(idxState,0,0,1))
    ]/terminalMass;

    terminalMeanMC = monte_carlo_terminal_mean( ...
        X0validation,T,sigma,rho,beta,center,scale,optsODE);

    fprintf('  terminal mass = %.8f\n', terminalMass);
    fprintf('  terminal mean = [% .6f % .6f % .6f]\n', terminalMean);
    fprintf('  Monte Carlo   = [% .6f % .6f % .6f]\n', terminalMeanMC);
    fprintf('  mean error    = %.3e\n', norm(terminalMean-terminalMeanMC));
end

%% Reconstruct terminal densities on the cube and evaluate on attractor

densityOnAttractor = cell(numel(terminalTimes),1);

for kk = 1:numel(terminalTimes)

    T = terminalTimes(kk);

    fprintf('\nReconstructing density for T = %g\n', T);

    [densityCoeff,densityBasis] = reconstruct_density_cube_3d( ...
        terminalMoments{kk}, basisState, momentDegree);

    rhoAttr = zeros(size(Yattr,1),1);

    for j = 1:size(Yattr,1)
        rhoAttr(j) = eval_poly_coeff_3d( ...
            densityCoeff,densityBasis, ...
            Yattr(j,1),Yattr(j,2),Yattr(j,3));
    end

    rhoClip = max(rhoAttr,0);
    rhoFloor = max(prctile_no_toolbox(rhoClip(rhoClip > 0),1),1e-12);

    densityOnAttractor{kk} = max(rhoClip,rhoFloor);

    fprintf('  raw density min/max = %.4e / %.4e\n', ...
        min(rhoAttr), max(rhoAttr));
end

%% Plot attractor colored by reconstructed terminal densities

for kk = 1:numel(terminalTimes)

    T = terminalTimes(kk);

    fig = figure('Color','w','Position',[100 100 1100 900]);

    ax = axes('Parent',fig);
    hold(ax,'on');

    ax.Position = [0.08 0.13 0.68 0.78];

    colorData = densityOnAttractor{kk};

    surface(ax, ...
        [Xraw(:,1) Xraw(:,1)], ...
        [Xraw(:,2) Xraw(:,2)], ...
        [Xraw(:,3) Xraw(:,3)], ...
        [colorData colorData], ...
        'FaceColor','none', ...
        'EdgeColor','interp', ...
        'LineWidth',5);

    %% Plot initial ball in original coordinates

    Y0plot = sample_ball(numInitialPlotPoints, ...
        initialCenterScaled(:).', initialRadiusScaled);

    X0plot = Y0plot.*scale + center;

    scatter3(ax, ...
        X0plot(:,1), X0plot(:,2), X0plot(:,3), ...
        18, 'k', 'filled', ...
        'MarkerFaceAlpha',1);

    colormap(ax,parula);

    cb = colorbar(ax);
    cb.Position = [0.84 0.18 0.035 0.68];
    cb.FontSize = 32;

    view(ax,35,22);

    set(ax,'FontSize',46);

    xlabel(ax,'$x$','Interpreter','latex','FontSize',56,'FontWeight','bold');
    ylabel(ax,'$y$','Interpreter','latex','FontSize',56,'FontWeight','bold');
    zlabel(ax,'$z$','Interpreter','latex','FontSize',56,'FontWeight','bold');

    axis(ax,'equal');
    axis(ax,'tight');
    grid(ax,'on');
    box(ax,'on');

    if saveFigures
        fileName = sprintf('lorenz_occ_density_T_%s.pdf', time_to_filename(T));
        save_pdf_figure(fig,fileName);
    end
end

fprintf('\nFinished Lorenz occupation-measure transport example.\n');

%% Local functions

function dx = lorenz_rhs(~,x,sigma,rho,beta)

    dx = [
        sigma*(x(2)-x(1));
        x(1)*(rho-x(3))-x(2);
        x(1)*x(2)-beta*x(3)
    ];
end

function F = scaled_lorenz_polynomial_terms(center,scale,sigma,rho,beta)
% Scaled vector field in variables Y = (X,Y,Z), where
%
%     original = center + scale .* Y.
%
% Each row has the form
%
%     [coefficient, exponent_X, exponent_Y, exponent_Z].

    c = center;
    s = scale;

    F = cell(3,1);

    F{1} = [
        sigma*(c(2)-c(1))/s(1), 0, 0, 0;
        -sigma,                 1, 0, 0;
        sigma*s(2)/s(1),        0, 1, 0
    ];

    F{2} = [
        (c(1)*(rho-c(3))-c(2))/s(2), 0, 0, 0;
        s(1)*(rho-c(3))/s(2),        1, 0, 0;
        -c(1)*s(3)/s(2),             0, 0, 1;
        -1,                          0, 1, 0;
        -s(1)*s(3)/s(2),             1, 0, 1
    ];

    F{3} = [
        (c(1)*c(2)-beta*c(3))/s(3), 0, 0, 0;
        s(1)*c(2)/s(3),             1, 0, 0;
        c(1)*s(2)/s(3),             0, 1, 0;
        s(1)*s(2)/s(3),             1, 1, 0;
        -beta,                      0, 0, 1
    ];
end

function constraints = build_liouville_constraints( ...
    yOcc, yTerm, yInitial, basisOcc, idxOcc, idxState, Fpoly, d, T)

    degf = 2;
    testDegree = d + 1 - degf;

    basisTest = basisOcc(sum(basisOcc,2) <= testDegree,:);

    constraints = [];

    for k = 1:size(basisTest,1)

        exponent = basisTest(k,:);
        kt = exponent(1);
        alpha = exponent(2:4);

        expr = 0;

        if kt > 0
            expr = expr + kt*moment_expr_4d( ...
                yOcc,idxOcc,kt-1,alpha(1),alpha(2),alpha(3));
        end

        if alpha(1) > 0
            expr = expr + T*alpha(1)*apply_vector_field_terms_4d( ...
                yOcc,idxOcc,kt,alpha-[1 0 0],Fpoly{1});
        end

        if alpha(2) > 0
            expr = expr + T*alpha(2)*apply_vector_field_terms_4d( ...
                yOcc,idxOcc,kt,alpha-[0 1 0],Fpoly{2});
        end

        if alpha(3) > 0
            expr = expr + T*alpha(3)*apply_vector_field_terms_4d( ...
                yOcc,idxOcc,kt,alpha-[0 0 1],Fpoly{3});
        end

        rhs = yTerm(get_idx_3d( ...
            idxState,alpha(1),alpha(2),alpha(3)));

        if kt == 0
            rhs = rhs - yInitial(get_idx_3d( ...
                idxState,alpha(1),alpha(2),alpha(3)));
        end

        constraints = [constraints, expr == rhs]; %#ok<AGROW>
    end
end

function expr = apply_vector_field_terms_4d(y,idx,kt,baseExp,terms)

    expr = 0;

    for r = 1:size(terms,1)

        coeff = terms(r,1);
        expv = baseExp + terms(r,2:4);

        expr = expr + coeff*moment_expr_4d( ...
            y,idx,kt,expv(1),expv(2),expv(3));
    end
end

function basis = monomial_basis_3d(d)

    basis = [];

    for total = 0:d
        for a = 0:total
            for b = 0:(total-a)
                c = total - a - b;
                basis = [basis; a b c]; %#ok<AGROW>
            end
        end
    end
end

function basis = monomial_basis_4d(d)

    basis = [];

    for total = 0:d
        for kt = 0:total
            for a = 0:(total-kt)
                for b = 0:(total-kt-a)
                    c = total - kt - a - b;
                    basis = [basis; kt a b c]; %#ok<AGROW>
                end
            end
        end
    end
end

function idx = make_index_map_3d(basis)

    idx = containers.Map();

    for k = 1:size(basis,1)
        idx(sprintf('%d_%d_%d',basis(k,1),basis(k,2),basis(k,3))) = k;
    end
end

function idx = make_index_map_4d(basis)

    idx = containers.Map();

    for k = 1:size(basis,1)
        idx(sprintf('%d_%d_%d_%d',basis(k,1),basis(k,2), ...
            basis(k,3),basis(k,4))) = k;
    end
end

function k = get_idx_3d(idx,a,b,c)

    key = sprintf('%d_%d_%d',a,b,c);

    if ~isKey(idx,key)
        error('Moment (%d,%d,%d) exceeds truncation degree.',a,b,c);
    end

    k = idx(key);
end

function k = get_idx_4d(idx,kt,a,b,c)

    key = sprintf('%d_%d_%d_%d',kt,a,b,c);

    if ~isKey(idx,key)
        error('Moment (%d,%d,%d,%d) exceeds truncation degree.', ...
            kt,a,b,c);
    end

    k = idx(key);
end

function expr = moment_expr_3d(y,idx,a,b,c)

    expr = y(get_idx_3d(idx,a,b,c));
end

function expr = moment_expr_4d(y,idx,kt,a,b,c)

    expr = y(get_idx_4d(idx,kt,a,b,c));
end

function M = build_moment_matrix_3d(y,idx,basisMoment)

    n = size(basisMoment,1);
    M = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisMoment(i,:) + basisMoment(j,:);

            M(i,j) = moment_expr_3d(y,idx,a(1),a(2),a(3));
            M(j,i) = M(i,j);
        end
    end
end

function M = build_moment_matrix_4d(y,idx,basisMoment)

    n = size(basisMoment,1);
    M = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisMoment(i,:) + basisMoment(j,:);

            M(i,j) = moment_expr_4d( ...
                y,idx,a(1),a(2),a(3),a(4));
            M(j,i) = M(i,j);
        end
    end
end

function Mg = build_cube_localizing_matrix_3d(y,idx,basisLocalizing,coord)

    n = size(basisLocalizing,1);
    Mg = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisLocalizing(i,:) + basisLocalizing(j,:);
            ashift = a;
            ashift(coord) = ashift(coord) + 2;

            val = moment_expr_3d(y,idx,a(1),a(2),a(3)) ...
                - moment_expr_3d(y,idx,ashift(1),ashift(2),ashift(3));

            Mg(i,j) = val;
            Mg(j,i) = val;
        end
    end
end

function Mg = build_cube_localizing_matrix_4d( ...
    y,idx,basisLocalizing,coord)

    n = size(basisLocalizing,1);
    Mg = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisLocalizing(i,:) + basisLocalizing(j,:);
            ashift = a;
            ashift(coord+1) = ashift(coord+1) + 2;

            val = moment_expr_4d(y,idx,a(1),a(2),a(3),a(4)) ...
                - moment_expr_4d(y,idx,ashift(1),ashift(2), ...
                    ashift(3),ashift(4));

            Mg(i,j) = val;
            Mg(j,i) = val;
        end
    end
end

function Mt = build_time_localizing_matrix_4d( ...
    y,idx,basisLocalizing)

    n = size(basisLocalizing,1);
    Mt = sdpvar(n,n,'symmetric');

    for i = 1:n
        for j = i:n

            a = basisLocalizing(i,:) + basisLocalizing(j,:);
            a1 = a;
            a2 = a;
            a1(1) = a1(1) + 1;
            a2(1) = a2(1) + 2;

            val = moment_expr_4d( ...
                y,idx,a1(1),a1(2),a1(3),a1(4)) ...
                - moment_expr_4d( ...
                y,idx,a2(1),a2(2),a2(3),a2(4));

            Mt(i,j) = val;
            Mt(j,i) = val;
        end
    end
end

function [coeff,basis] = reconstruct_density_cube_3d(yval,basisY,d)

    basis = monomial_basis_3d(d);
    N = size(basis,1);

    idxY = make_index_map_3d(basisY);

    massMatrix = zeros(N,N);
    rhs = zeros(N,1);

    for i = 1:N

        a = basis(i,:);

        rhs(i) = yval(get_idx_3d(idxY,a(1),a(2),a(3)));

        for j = 1:N
            p = basis(i,:) + basis(j,:);
            massMatrix(i,j) = cube_moment_3d(p(1),p(2),p(3));
        end
    end

    reg = 1e-10;
    coeff = (massMatrix + reg*eye(N))\rhs;
end

function val = cube_moment_3d(a,b,c)

    val = one_dim_cube_moment(a) ...
        * one_dim_cube_moment(b) ...
        * one_dim_cube_moment(c);
end

function val = one_dim_cube_moment(k)

    if mod(k,2) == 1
        val = 0;
    else
        val = 2/(k+1);
    end
end

function val = eval_poly_coeff_3d(coeff,basis,x,y,z)

    val = 0;

    for k = 1:length(coeff)
        val = val + coeff(k)*x^basis(k,1)*y^basis(k,2)*z^basis(k,3);
    end
end

function m = shifted_ball_moment(alpha,center,r)

    m = 0;

    for k1 = 0:alpha(1)
        for k2 = 0:alpha(2)
            for k3 = 0:alpha(3)

                k = [k1 k2 k3];

                binom = nchoosek(alpha(1),k1) ...
                    * nchoosek(alpha(2),k2) ...
                    * nchoosek(alpha(3),k3);

                centerPart = prod(center.^(alpha-k));

                m = m + binom*centerPart*central_ball_moment(k,r);
            end
        end
    end
end

function m = central_ball_moment(k,r)
% Normalized moments of the uniform probability measure on a centered
% radius-r ball in R^3.

    if any(mod(k,2) == 1)
        m = 0;
        return;
    end

    p = k/2;
    K = sum(p);

    m = r^(2*K)*gamma(3/2 + 1)/gamma(K + 3/2 + 1);

    for i = 1:3
        m = m*gamma(p(i)+1/2)/gamma(1/2);
    end
end

function X = sample_ball(N,center,r)

    U = randn(N,3);
    U = U./vecnorm(U,2,2);

    R = r*rand(N,1).^(1/3);

    X = center + U.*R;
end

function q = prctile_no_toolbox(x,p)

    x = sort(x(:));

    if isempty(x)
        q = 0;
        return;
    end

    k = 1 + (numel(x)-1)*p/100;
    k0 = floor(k);
    k1 = ceil(k);

    if k0 == k1
        q = x(k0);
    else
        q = x(k0) + (k-k0)*(x(k1)-x(k0));
    end
end

function meanScaled = monte_carlo_terminal_mean( ...
    X0,T,sigma,rho,beta,center,scale,optsODE)

    n = size(X0,1);
    finalScaled = zeros(n,3);

    for j = 1:n
        [~,trajectory] = ode45( ...
            @(t,x) lorenz_rhs(t,x,sigma,rho,beta), ...
            [0 T],X0(j,:).',optsODE);
        finalScaled(j,:) = (trajectory(end,:)-center)./scale;
    end

    meanScaled = mean(finalScaled,1).';
end

function str = time_to_filename(T)

    str = sprintf('%.10g',T);
    str = strrep(str,'.','p');
    str = strrep(str,'-','m');
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