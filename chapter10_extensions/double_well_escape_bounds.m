%% double_well_escape_bounds.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Extensions and Frontiers
% Example: Escape rates for Brownian dynamics in a double-well potential
%
% This script computes SOS bounds on the principal eigenvalue associated
% with escape from a one-dimensional double-well potential.
%
% The computation is performed on the rescaled interval
%
%     -1 <= x <= 1,
%
% where the physical coordinate is X = L x.
%
% Pressing Run should reproduce the numerical output and save:
%
%     double_well_potential.pdf
%     double_well_eigenfunction.pdf
%
% Requirements:
%   - MATLAB
%   - YALMIP
%   - MOSEK
%
% -------------------------------------------------------------------------

clear; clc; close all;
yalmip('clear');
format long;

%% User options

sigmaNoise = 1.0;
alpha = 0.75;
L = 3.0;

degree = 10;

xEnd = 1 - 1e-4;
numGrid = 2001;

saveFigures = true;

%% Solver options

opts = sdpsettings( ...
    'solver','mosek', ...
    'verbose',1, ...
    'sos.model',2);

%% Book plotting style

S.black = [0, 0, 0];
S.blue  = [0, 92, 175]/255;

%% Symbolic variables

sdpvar lambda x u v

%% Witten-Laplacian potential in rescaled coordinates

w = -0.5*( ...
        -alpha*(2*L*x - 1) ...
        + (2 - 3*L*x)*L*x ) ...
    + (1/(sigmaNoise^2))*((L*x)*(1 - L*x)*(L*x + alpha))^2;

%% Variational SOS inequality

p = (sigmaNoise^2)*v^2/L^2 + w*u^2 - lambda*u^2;

[f,cf] = polynomial(x,degree);

dfdx = jacobian(f,x);
df = dfdx*(1 - x^2) + 2*f*x;

%% Polynomial multiplier with degree restriction in u and v

powers = monpowers(3,degree+2);
powers = powers(powers(:,2) + powers(:,3) <= 2,:);

cs = sdpvar(size(powers,1),1);

s = 0;

for k = 1:size(powers,1)
    s = s + cs(k)*prod([x u v].^powers(k,:));
end

%% SOS lower-bound problem

sosPolynomial = ...
    p*(1 - x^2)^2 ...
    + df*u^2 ...
    + 2*f*v*u*(1 - x^2) ...
    - (1 - x^2)*s;

constraints = [
    sos(sosPolynomial), ...
    sos(s)
    ];

fprintf('\n============================================================\n');
fprintf('Double-well escape-rate SOS bound\n');
fprintf('============================================================\n\n');

fprintf('sigma = %.6g\n',sigmaNoise);
fprintf('alpha = %.6g\n',alpha);
fprintf('L     = %.6g\n',L);
fprintf('degree = %d\n\n',degree);

fprintf('Solving SOS eigenvalue lower-bound problem...\n');

sol = solvesos(constraints,-lambda,opts,[cf; cs; lambda]);

if sol.problem ~= 0
    warning('SOS solve returned status: %s',sol.info);
end

lambdaLower = value(lambda);

fprintf('SOS lower bound on principal eigenvalue: %.12e\n',lambdaLower);

%% Recover approximate eigenfunction

xPlot = linspace(-xEnd,xEnd,numGrid);

fApprox = clean(replace(f,cf,value(cf)),1e-10);
fHandle = yalmip_polynomial_to_function(fApprox,x);

t = linspace(0,1,2001);

logU = zeros(size(xPlot));

for i = 1:numel(xPlot)

    xi = xPlot(i);
    xt = xi*t;

    logGradU = xi*fHandle(xt)./(1 - xt.^2);

    logU(i) = -trapz(t,logGradU);
end

uStar = exp(logU*L^2/sigmaNoise^2);

%% Undo symmetrization

physicalX = L*xPlot;

potential = double_well_potential(physicalX,alpha);

vStar = uStar.*exp(-0.5*potential/(sigmaNoise^2/L^2));

%% Rayleigh quotient upper bound

normU = trapz(xPlot,uStar.^2);

gradU = uStar.*fHandle(xPlot)./(1 - xPlot.^2);

normGradU = (L^2/sigmaNoise^2)*trapz(xPlot,gradU.^2);

wHandle = yalmip_polynomial_to_function(w,x);

normPotential = trapz(xPlot,wHandle(xPlot).*uStar.^2);

lambdaUpper = (normGradU + normPotential)/normU;

%% Print summary

fprintf('\nBounds using degree %d polynomials:\n',degree);
fprintf('  Upper bound = %.12e\n',lambdaUpper);
fprintf('  Lower bound = %.12e\n',lambdaLower);
fprintf('  Gap         = %.12e\n',lambdaUpper - lambdaLower);

%% Figure 1: potential landscape

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 11 8]);

plot(physicalX,potential, ...
    'Color',S.black, ...
    'LineWidth',3);

xlabel('$x$','Interpreter','latex','FontSize',28,'FontWeight','bold');
ylabel('$V(x)$','Interpreter','latex','FontSize',28,'FontWeight','bold');

set(gca,'FontSize',22,'TickLabelInterpreter','latex');

box on
grid on

if saveFigures
    export_pdf(fig1,'double_well_potential.pdf');
end

%% Figure 2: approximate leading eigenfunction

fig2 = figure;
set(fig2,'Color','w');

plot(physicalX,vStar, ...
    'Color',S.blue, ...
    'LineWidth',4);

xlabel('$x$','Interpreter','latex','FontSize',36,'FontWeight','bold');
ylabel('$u_0(x)$','Interpreter','latex','FontSize',36,'FontWeight','bold');

set(gca,'FontSize',24,'TickLabelInterpreter','latex');

box on
grid on
xlim([-3 3])

if saveFigures
    export_pdf(fig2,'double_well_eigenfunction.pdf');
end

fprintf('\nFinished double-well escape-rate example.\n');

%% Local functions

function V = double_well_potential(x,alpha)

    V = (1/12)*x.^2.*( ...
        -alpha*(6 - 4*x) ...
        + x.*(3*x - 4));
end

function fh = yalmip_polynomial_to_function(poly,var)

    str = sdisplay(poly);
    str = str{1};

    varName = sdisplay(var);
    varName = varName{1};

    str = strrep(str,'^','.^');
    str = strrep(str,'*','.*');

    fh = str2func(['@(' varName ')' str]);
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