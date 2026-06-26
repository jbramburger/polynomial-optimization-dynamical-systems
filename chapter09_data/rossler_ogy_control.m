%% rossler_ogy_control.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Data-Driven Auxiliary Function Methods
% Example: Data-driven OGY control of Rössler periodic orbits
%
% This script stabilizes unstable periodic orbits of the Rössler Poincaré
% map using data-driven OGY control.
%
% The Rössler parameter c is used as a scalar control input:
%
%     c = c0 + u.
%
% The script:
%
%   1. collects Poincaré map data for c near c0 = 18,
%   2. identifies a fixed point and a period-two orbit from uncontrolled data,
%   3. learns local control-affine return maps from data,
%   4. computes OGY pole-placement gains,
%   5. certifies local Lyapunov activation intervals by grid search, and
%   6. simulates the controlled full Rössler flow.
%
% Pressing Run should reproduce the numerical output and save:
%
%   rossler_control_period_1.pdf
%   rossler_control_period_2.pdf
%   rossler_ogy_control_results.mat
%
% Requirements:
%   - MATLAB
%
% -------------------------------------------------------------------------

clear; clc; close all;
rng(1);

%% User options

a  = 0.1;
b  = 0.1;
c0 = 18.0;

cVals = c0 + linspace(-0.25,0.25,11);

numTransientCrossings = 150;
numDataCrossings = 900;

uMax = 0.20;

degMap = 7;
mapReg = 1e-10;

lambdaDesiredFixed = 0.20;
lambdaDesiredPeriodTwo = 0.20;

fixedDataRadius = 2.50;
periodDataRadius = 1.25;

certTol = 1e-9;
qSearchMaxFixed = 1.25;
qSearchMaxPeriodTwo = 1.25;
numCertGrid = 5001;
numDeltaGrid = 2000;

numControlCrossingsFixed = 300;
numControlCrossingsPeriodTwo = 400;

saveFigures = true;
saveData = true;

%% Plot style

S.blue   = [0, 92, 175]/255;
S.red    = [200, 50, 50]/255;
S.black  = [0, 0, 0];

%% ODE options

odeOpts = odeset( ...
    'RelTol',1e-10, ...
    'AbsTol',1e-12, ...
    'MaxStep',0.05, ...
    'Events',@(t,x) rossler_section_event(t,x));

%% Collect controlled Poincare map data

fprintf('\n============================================================\n');
fprintf('Data-driven OGY control for the Rössler Poincaré map\n');
fprintf('============================================================\n\n');

fprintf('Collecting Poincare map data...\n');

Y0All = [];
Y1All = [];
UAll  = [];
CAll  = [];

crossingStatesByC = cell(numel(cVals),1);

for r = 1:numel(cVals)

    c = cVals(r);

    fprintf('  c = %.4f\n',c);

    xInit = [1; 0; 0];

    states = collect_rossler_crossings( ...
        xInit,a,b,c, ...
        numTransientCrossings + numDataCrossings + 1, ...
        odeOpts);

    states = states(numTransientCrossings+1:end,:);
    crossingStatesByC{r} = states;

    y0 = states(1:end-1,2);
    y1 = states(2:end,2);

    Y0All = [Y0All; y0]; %#ok<AGROW>
    Y1All = [Y1All; y1]; %#ok<AGROW>
    UAll  = [UAll; (c-c0)*ones(size(y0))]; %#ok<AGROW>
    CAll  = [CAll; c*ones(size(y0))]; %#ok<AGROW>
end

numTriples = numel(Y0All);

fprintf('Total snapshot triples: %d\n',numTriples);

%% Uncontrolled return map data

idxUncontrolled = abs(CAll-c0) < 1e-12;

y0Uncontrolled = Y0All(idxUncontrolled);
y1Uncontrolled = Y1All(idxUncontrolled);

[~,ysu,fyu] = build_scalar_interpolant(y0Uncontrolled,y1Uncontrolled);

%% Fixed point identification

ystar = find_fixed_point_from_data(ysu,fyu);

fprintf('\n============================================================\n');
fprintf('Fixed-point OGY stabilization\n');
fprintf('============================================================\n');
fprintf('Estimated fixed point y_* = %.12f\n',ystar);

Q0 = Y0All - ystar;
Q1 = Y1All - ystar;

idxFixedLocal = abs(Q0) <= fixedDataRadius;

[F0Fixed,F1Fixed,fitErrFixed] = learn_controlled_scalar_map( ...
    Q0(idxFixedLocal), ...
    Q1(idxFixedLocal), ...
    UAll(idxFixedLocal), ...
    degMap,mapReg,1.25,true);

fprintf('Fixed-point controlled map fit residual = %.4e\n',fitErrFixed);

AFixed = F0Fixed(2);
BFixed = F1Fixed(1);

if abs(BFixed) < 1e-10
    error('Fixed-point control sensitivity is too small.');
end

KFixed = (lambdaDesiredFixed - AFixed)/BFixed;

fprintf('Local learned linearization:\n');
fprintf('  A       = %.8e\n',AFixed);
fprintf('  B       = %.8e\n',BFixed);
fprintf('  K       = %.8e\n',KFixed);
fprintf('  A + BK  = %.8e\n',AFixed + BFixed*KFixed);

[deltaCertFixed,certInfoFixed] = maximize_lyapunov_interval( ...
    F0Fixed,F1Fixed,KFixed,uMax, ...
    qSearchMaxFixed,numCertGrid,numDeltaGrid,certTol);

deltaSatFixed = 0.9*uMax/max(abs(KFixed),1e-12);
qActFixed = min(deltaCertFixed,deltaSatFixed);

if qActFixed <= 0
    warning('No certified fixed-point interval found. Falling back to small activation radius.');
    qActFixed = min(0.05,deltaSatFixed);
end

fprintf('Certified Lyapunov interval delta = %.8e\n',deltaCertFixed);
fprintf('Saturation-limited interval       = %.8e\n',deltaSatFixed);
fprintf('Final activation radius qAct      = %.8e\n',qActFixed);

%% Simulate fixed-point controlled system

fprintf('\nSimulating fixed-point controlled Rössler system...\n');

[~,idxClosest] = min(abs(cVals-c0));
statesUncontrolled = crossingStatesByC{idxClosest};

xCross = statesUncontrolled(10,:).';

[tContinuousFixed,yContinuousFixed,tSectionFixed,ySectionFixed, ...
    YclFixed,QclFixed,UclFixed,CclFixed,VclFixed,activeFixed] = ...
    simulate_fixed_point_control( ...
        xCross,a,b,c0,ystar,KFixed,qActFixed,uMax, ...
        numControlCrossingsFixed);

fprintf('Fixed-point control activated on %d / %d crossings.\n', ...
    sum(activeFixed),numControlCrossingsFixed);

%% Book figure: fixed point

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 13 7]);

plot(tContinuousFixed,yContinuousFixed,'k','LineWidth',3);
hold on
plot(tSectionFixed,ySectionFixed,'.', ...
    'Color',S.red, ...
    'MarkerSize',30, ...
    'LineWidth',1.2);

set(gca,'FontSize',18,'TickLabelInterpreter','latex');

xlabel('$t$','Interpreter','latex','FontSize',24);
ylabel('$y(t)$','Interpreter','latex','FontSize',24);

xlim([0 200]);
ylim([-29 25]);

box on
grid on

if saveFigures
    export_pdf(fig1,'rossler_control_period_1.pdf');
end

%% Period-two orbit identification

fprintf('\n============================================================\n');
fprintf('Period-two OGY stabilization\n');
fprintf('============================================================\n');

[y2star,periodTwoFound] = find_period_two_orbit_from_data( ...
    y0Uncontrolled,y1Uncontrolled,ystar);

if ~periodTwoFound

    warning('Could not reliably identify a period-two orbit from the data.');

    ctrlPeriodTwo = [];
    tContinuousPeriodTwo = [];
    yContinuousPeriodTwo = [];
    tSectionPeriodTwo = [];
    ySectionPeriodTwo = [];

else

    fprintf('Estimated period-two points:\n');
    fprintf('  y_1^* = %.12f\n',y2star(1));
    fprintf('  y_2^* = %.12f\n',y2star(2));

    ctrlPeriodTwo(1) = build_phase_controller( ...
        Y0All,Y1All,UAll, ...
        y2star(1),y2star(2), ...
        degMap,mapReg,periodDataRadius, ...
        lambdaDesiredPeriodTwo,uMax, ...
        qSearchMaxPeriodTwo,numCertGrid,numDeltaGrid,certTol);

    ctrlPeriodTwo(2) = build_phase_controller( ...
        Y0All,Y1All,UAll, ...
        y2star(2),y2star(1), ...
        degMap,mapReg,periodDataRadius, ...
        lambdaDesiredPeriodTwo,uMax, ...
        qSearchMaxPeriodTwo,numCertGrid,numDeltaGrid,certTol);

    for j = 1:2

        fprintf('\nPhase %d controller:\n',j);
        fprintf('  yFrom = %.12f\n',ctrlPeriodTwo(j).yFrom);
        fprintf('  yTo   = %.12f\n',ctrlPeriodTwo(j).yTo);
        fprintf('  fit residual = %.4e\n',ctrlPeriodTwo(j).fitErr);
        fprintf('  A = %.8e\n',ctrlPeriodTwo(j).A);
        fprintf('  B = %.8e\n',ctrlPeriodTwo(j).B);
        fprintf('  K = %.8e\n',ctrlPeriodTwo(j).K);
        fprintf('  A+B*K = %.8e\n',ctrlPeriodTwo(j).A + ctrlPeriodTwo(j).B*ctrlPeriodTwo(j).K);
        fprintf('  certified delta = %.8e\n',ctrlPeriodTwo(j).deltaCert);
        fprintf('  qAct = %.8e\n',ctrlPeriodTwo(j).qAct);
    end

    %% Simulate period-two controlled system

    fprintf('\nSimulating period-two controlled Rössler system...\n');

    xCross = statesUncontrolled(20,:).';

    [tContinuousPeriodTwo,yContinuousPeriodTwo,tSectionPeriodTwo,ySectionPeriodTwo, ...
        YclPeriodTwo,UclPeriodTwo,CclPeriodTwo,phasePeriodTwo, ...
        distPeriodTwo,activePeriodTwo] = ...
        simulate_period_two_control( ...
            xCross,a,b,c0,ctrlPeriodTwo,uMax, ...
            numControlCrossingsPeriodTwo);

    fprintf('Period-two control activated on %d / %d crossings.\n', ...
        sum(activePeriodTwo),numControlCrossingsPeriodTwo);

    %% Book figure: period two

    fig2 = figure;
    set(fig2,'Color','w','Units','centimeters','Position',[2 2 13 7]);

    plot(tContinuousPeriodTwo,yContinuousPeriodTwo,'k','LineWidth',3);
    hold on
    plot(tSectionPeriodTwo,ySectionPeriodTwo,'.', ...
        'Color',S.blue, ...
        'MarkerSize',30, ...
        'LineWidth',1.2);

    set(gca,'FontSize',18,'TickLabelInterpreter','latex');

    xlabel('$t$','Interpreter','latex','FontSize',24);
    ylabel('$y(t)$','Interpreter','latex','FontSize',24);

    xlim([0 200]);
    ylim([-29 25]);

    box on
    grid on

    if saveFigures
        export_pdf(fig2,'rossler_control_period_2.pdf');
    end
end

%% Save data

if saveData

    save('rossler_ogy_control_results.mat', ...
        'a','b','c0','cVals','Y0All','Y1All','UAll','CAll', ...
        'ystar','F0Fixed','F1Fixed','AFixed','BFixed','KFixed', ...
        'deltaCertFixed','deltaSatFixed','qActFixed','certInfoFixed', ...
        'tContinuousFixed','yContinuousFixed','tSectionFixed','ySectionFixed', ...
        'YclFixed','QclFixed','UclFixed','CclFixed','VclFixed','activeFixed', ...
        'y2star','periodTwoFound','ctrlPeriodTwo', ...
        'tContinuousPeriodTwo','yContinuousPeriodTwo', ...
        'tSectionPeriodTwo','ySectionPeriodTwo');

    fprintf('\nSaved data: rossler_ogy_control_results.mat\n');
end

fprintf('\nFinished Rössler OGY control example.\n');

%% Local functions

function dx = rossler_rhs(~,x,a,b,c)

    dx = zeros(3,1);

    dx(1) = -x(2) - x(3);
    dx(2) =  x(1) + a*x(2);
    dx(3) =  b + x(3)*(x(1)-c);
end

function [value,isterminal,direction] = rossler_section_event(~,x)

    value = x(1);
    isterminal = 1;
    direction = +1;
end

function states = collect_rossler_crossings(x0,a,b,c,numCrossings,odeOpts)

    states = zeros(numCrossings,3);

    x = x0(:);

    [~,xx] = ode45(@(t,x) rossler_rhs(t,x,a,b,c),[0 200],x, ...
        odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.05));

    x = xx(end,:).';

    for k = 1:numCrossings
        x = next_rossler_crossing(x,a,b,c,odeOpts);
        states(k,:) = x.';
    end
end

function xNext = next_rossler_crossing(x0,a,b,c,odeOpts)

    tinyOpts = odeset( ...
        'RelTol',1e-10, ...
        'AbsTol',1e-12, ...
        'MaxStep',1e-3);

    [~,xx] = ode45(@(t,x) rossler_rhs(t,x,a,b,c),[0 1e-3],x0(:),tinyOpts);

    xStart = xx(end,:).';

    [~,~,te,xe] = ode45( ...
        @(t,x) rossler_rhs(t,x,a,b,c), ...
        [0 300], ...
        xStart, ...
        odeOpts);

    if isempty(te)
        error('No Poincare crossing detected. Increase tMax or check dynamics.');
    end

    xNext = xe(end,:).';
end

function [F,ysu,fyu] = build_scalar_interpolant(y0,y1)

    [ys,ord] = sort(y0(:));
    fy = y1(ord);

    [ysu,ia] = unique(ys,'stable');
    fyu = fy(ia);

    F = @(yy) interp1(ysu,fyu,yy,'pchip');
end

function ystar = find_fixed_point_from_data(ysu,fyu)

    fsu = fyu - ysu;

    rootCandidates = [];

    for j = 1:numel(ysu)-1

        if fsu(j) == 0

            rootCandidates(end+1) = ysu(j); %#ok<AGROW>

        elseif fsu(j)*fsu(j+1) < 0

            try
                rootCandidates(end+1) = fzero( ...
                    @(yy) interp1(ysu,fsu,yy,'pchip'), ...
                    [ysu(j),ysu(j+1)]); %#ok<AGROW>
            catch
            end
        end
    end

    if isempty(rootCandidates)

        [~,jmin] = min(abs(fsu));
        ystar = ysu(jmin);

        warning('No sign-change fixed point found. Using nearest data point.');

    else

        ystar = rootCandidates(1);
    end
end

function [F0,F1,relErr] = learn_controlled_scalar_map( ...
    q0,q1,u,degMap,lambdaReg,weightRadius,forceZero)

    q0 = q0(:);
    q1 = q1(:);
    u  = u(:);

    M = degMap + 1;

    X = zeros(M,numel(q0));

    for k = 1:numel(q0)
        X(:,k) = monomial_vector(q0(k),degMap);
    end

    A = [X; X.*u.'];

    w = exp(-(abs(q0)/weightRadius).^8);
    Wsqrt = sqrt(w(:)).';

    AW = A.*Wsqrt;
    YW = q1(:).'.*Wsqrt;

    Fall = YW*AW'/(AW*AW' + lambdaReg*eye(2*M));

    F0 = Fall(1:M);
    F1 = Fall(M+1:end);

    if forceZero
        F0(1) = 0;
    end

    relErr = norm(YW - Fall*AW,'fro')/max(1e-14,norm(YW,'fro'));
end

function [deltaCert,info] = maximize_lyapunov_interval( ...
    F0,F1,K,uMax,qSearchMax,numCertGrid,numDeltaGrid,certTol)

    qGrid = linspace(-qSearchMax,qSearchMax,numCertGrid).';
    qGrid = qGrid(abs(qGrid) > 1e-8);

    uGrid = saturate(K*qGrid,uMax);

    qNext = eval_poly_row(F0,qGrid) + eval_poly_row(F1,qGrid).*uGrid;

    deltaV = qNext.^2 - qGrid.^2;

    deltaCandidates = linspace(0.001,qSearchMax,numDeltaGrid);

    deltaCert = 0;

    for j = 1:numel(deltaCandidates)

        delta = deltaCandidates(j);
        idx = abs(qGrid) <= delta;

        if all(deltaV(idx) <= certTol)
            deltaCert = delta;
        else
            break;
        end
    end

    info = struct();
    info.qGrid = qGrid;
    info.uGrid = uGrid;
    info.qNext = qNext;
    info.deltaV = deltaV;
end

function [tCont,yCont,tSec,ySec,Ycl,Qcl,Ucl,Ccl,Vcl,active] = ...
    simulate_fixed_point_control( ...
        xCross,a,b,c0,ystar,K,qAct,uMax,numCrossings)

    tNow = 0;

    tCont = [];
    yCont = [];

    tSec = zeros(numCrossings,1);
    ySec = zeros(numCrossings,1);

    Ycl = zeros(numCrossings,1);
    Qcl = zeros(numCrossings,1);
    Ucl = zeros(numCrossings,1);
    Ccl = zeros(numCrossings,1);
    Vcl = zeros(numCrossings,1);
    active = false(numCrossings,1);

    for k = 1:numCrossings

        yNow = xCross(2);
        qNow = yNow - ystar;

        if abs(qNow) <= qAct
            uNow = saturate(K*qNow,uMax);
            active(k) = true;
        else
            uNow = 0;
            active(k) = false;
        end

        cNow = c0 + uNow;

        Ycl(k) = yNow;
        Qcl(k) = qNow;
        Ucl(k) = uNow;
        Ccl(k) = cNow;
        Vcl(k) = qNow^2;

        [xCross,tSeg,xSeg,tEvent,xEvent] = ...
            next_rossler_crossing_with_trajectory(xCross,a,b,cNow,tNow);

        if isempty(tCont)
            tCont = tSeg;
            yCont = xSeg(:,2);
        else
            tCont = [tCont; tSeg(2:end)]; %#ok<AGROW>
            yCont = [yCont; xSeg(2:end,2)]; %#ok<AGROW>
        end

        tNow = tEvent;
        tSec(k) = tEvent;
        ySec(k) = xEvent(2);
    end
end

function [y2star,success] = find_period_two_orbit_from_data(y0,y1,ystar)

    success = false;
    y2star = [NaN; NaN];

    [F,ysu,~] = build_scalar_interpolant(y0,y1);

    ymin = min(ysu);
    ymax = max(ysu);

    grid = linspace(ymin,ymax,5000).';
    Fgrid = F(grid);

    valid1 = isfinite(Fgrid) & Fgrid >= ymin & Fgrid <= ymax;

    H = NaN(size(grid));
    H(valid1) = F(Fgrid(valid1)) - grid(valid1);

    roots = [];

    for j = 1:numel(grid)-1

        if ~isfinite(H(j)) || ~isfinite(H(j+1))
            continue;
        end

        if H(j) == 0

            roots(end+1,1) = grid(j); %#ok<AGROW>

        elseif H(j)*H(j+1) < 0

            try
                rr = fzero(@(yy) F(F(yy)) - yy,[grid(j),grid(j+1)]);
                roots(end+1,1) = rr; %#ok<AGROW>
            catch
            end
        end
    end

    if isempty(roots)
        return;
    end

    roots = sort(roots(:));

    keep = true(size(roots));

    for j = 2:numel(roots)
        if abs(roots(j)-roots(j-1)) < 1e-6
            keep(j) = false;
        end
    end

    roots = roots(keep);
    roots = roots(abs(roots-ystar) > 1e-2);

    if numel(roots) < 2
        return;
    end

    bestErr = inf;
    bestPair = [];

    for i = 1:numel(roots)

        yi = roots(i);
        image = F(yi);

        [err,j] = min(abs(roots-image));

        if j ~= i && err < bestErr
            bestErr = err;
            bestPair = [yi; roots(j)];
        end
    end

    if isempty(bestPair)
        return;
    end

    y2star = bestPair;

    if abs(F(y2star(1))-y2star(2)) > abs(F(y2star(2))-y2star(1))
        y2star = flipud(y2star);
    end

    success = true;
end

function ctrl = build_phase_controller( ...
    Y0,Y1,U,yFrom,yTo,degMap,mapReg,dataRadius, ...
    lambdaDes,uMax,qSearchMax,numCertGrid,numDeltaGrid,certTol)

    idx = abs(Y0-yFrom) <= dataRadius;

    q0 = Y0(idx) - yFrom;
    q1 = Y1(idx) - yTo;
    uu = U(idx);

    if numel(q0) < 20
        warning('Very few data points for phase controller.');
    end

    [F0,F1,fitErr] = learn_controlled_scalar_map( ...
        q0,q1,uu,degMap,mapReg,0.75*dataRadius,true);

    A = F0(2);
    B = F1(1);

    if abs(B) < 1e-10
        error('Local control sensitivity B is too small.');
    end

    K = (lambdaDes - A)/B;

    [deltaCert,certInfo] = maximize_lyapunov_interval( ...
        F0,F1,K,uMax,qSearchMax,numCertGrid,numDeltaGrid,certTol);

    deltaSat = 0.9*uMax/max(abs(K),1e-12);
    qAct = min(deltaCert,deltaSat);

    if qAct <= 0
        warning('No certified phase interval found. Falling back to small activation radius.');
        qAct = min(0.05,deltaSat);
    end

    ctrl = struct();
    ctrl.yFrom = yFrom;
    ctrl.yTo = yTo;
    ctrl.F0 = F0;
    ctrl.F1 = F1;
    ctrl.fitErr = fitErr;
    ctrl.A = A;
    ctrl.B = B;
    ctrl.K = K;
    ctrl.lambdaDes = lambdaDes;
    ctrl.deltaCert = deltaCert;
    ctrl.deltaSat = deltaSat;
    ctrl.qAct = qAct;
    ctrl.certInfo = certInfo;
end

function [u,phase] = period_two_ogy_control(y,ctrl2,uMax)

    d1 = abs(y - ctrl2(1).yFrom);
    d2 = abs(y - ctrl2(2).yFrom);

    if d1 <= d2
        phase = 1;
        dmin = d1;
    else
        phase = 2;
        dmin = d2;
    end

    if dmin <= ctrl2(phase).qAct

        q = y - ctrl2(phase).yFrom;
        u = saturate(ctrl2(phase).K*q,uMax);

    else

        u = 0;
        phase = 0;
    end
end

function [tCont,yCont,tSec,ySec,Ycl,Ucl,Ccl,phase,dist,active] = ...
    simulate_period_two_control( ...
        xCross,a,b,c0,ctrl2,uMax,numCrossings)

    tNow = 0;

    tCont = [];
    yCont = [];

    tSec = zeros(numCrossings,1);
    ySec = zeros(numCrossings,1);

    Ycl = zeros(numCrossings,1);
    Ucl = zeros(numCrossings,1);
    Ccl = zeros(numCrossings,1);
    phase = zeros(numCrossings,1);
    dist = zeros(numCrossings,1);
    active = false(numCrossings,1);

    y2star = [ctrl2(1).yFrom; ctrl2(2).yFrom];

    for k = 1:numCrossings

        yNow = xCross(2);

        [uNow,phaseNow] = period_two_ogy_control(yNow,ctrl2,uMax);

        cNow = c0 + uNow;

        Ycl(k) = yNow;
        Ucl(k) = uNow;
        Ccl(k) = cNow;
        phase(k) = phaseNow;
        active(k) = phaseNow > 0;
        dist(k) = min(abs(yNow - y2star(:)));

        [xCross,tSeg,xSeg,tEvent,xEvent] = ...
            next_rossler_crossing_with_trajectory(xCross,a,b,cNow,tNow);

        if isempty(tCont)
            tCont = tSeg;
            yCont = xSeg(:,2);
        else
            tCont = [tCont; tSeg(2:end)]; %#ok<AGROW>
            yCont = [yCont; xSeg(2:end,2)]; %#ok<AGROW>
        end

        tNow = tEvent;
        tSec(k) = tEvent;
        ySec(k) = xEvent(2);
    end
end

function v = monomial_vector(x,d)

    v = zeros(d+1,1);

    for j = 0:d
        v(j+1) = x.^j;
    end
end

function vals = eval_poly_row(coef,q)

    coef = coef(:);
    d = numel(coef)-1;

    vals = zeros(size(q));

    for j = 0:d
        vals = vals + coef(j+1)*q.^j;
    end
end

function y = saturate(x,xmax)

    y = max(min(x,xmax),-xmax);
end

function [xNext,tSeg,xSeg,tEvent,xEvent] = ...
    next_rossler_crossing_with_trajectory(x0,a,b,c,t0)

    tinyOpts = odeset( ...
        'RelTol',1e-10, ...
        'AbsTol',1e-12, ...
        'MaxStep',1e-3);

    [ttiny,xtiny] = ode45(@(t,x) rossler_rhs(t,x,a,b,c), ...
        [0 1e-3],x0(:),tinyOpts);

    xStart = xtiny(end,:).';

    eventOpts = odeset( ...
        'RelTol',1e-10, ...
        'AbsTol',1e-12, ...
        'MaxStep',0.05, ...
        'Events',@rossler_section_event);

    [t,x,te,xe] = ode45( ...
        @(t,x) rossler_rhs(t,x,a,b,c), ...
        [0 300], ...
        xStart, ...
        eventOpts);

    if isempty(te)
        error('No Poincare crossing detected. Increase tMax or check dynamics.');
    end

    xNext = xe(end,:).';
    tEvent = t0 + 1e-3 + te(end);
    xEvent = xe(end,:).';

    tSeg = [ttiny; 1e-3 + t];
    xSeg = [xtiny; x];

    tSeg = t0 + tSeg;
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