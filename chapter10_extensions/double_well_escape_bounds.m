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

%% Witten-Laplacian potential used by the variational problem

w = -0.5*( ...
        -alpha*(2*L*x - 1) ...
        + (2 - 3*L*x)*L*x ) ...
    + (1/(sigmaNoise^2))*((L*x)*(1 - L*x)*(L*x + alpha))^2;

%% Variational integrand and rational null-Lagrangian field

kappa = sigmaNoise^2/L^2;
p = kappa*v^2 + w*u^2 - lambda*u^2;

% The optimized rational field is f(x)=fNumerator(x)/(1-x^2).
[fNumerator,cf] = polynomial(x,degree);

dfdx = jacobian(fNumerator,x);
g = 1-x^2;

% g^2*d(fNumerator/g)/dx = g*fNumerator' + 2*x*fNumerator.
divergenceNumerator = dfdx*g + 2*fNumerator*x;

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
    p*g^2 ...
    + divergenceNumerator*u^2 ...
    + 2*fNumerator*v*u*g ...
    - g*s;

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

%% Reconstruct the eigenfunction from the optimized null Lagrangian

xPlot = linspace(-xEnd,xEnd,numGrid);

fNumeratorApprox = clean( ...
    replace(fNumerator,cf,value(cf)),1e-10);
fNumeratorHandle = yalmip_polynomial_to_function( ...
    fNumeratorApprox,x);

t = linspace(0,1,2001);
logEigenfunction = zeros(size(xPlot));

for i = 1:numel(xPlot)
    xi = xPlot(i);
    xt = xi*t;

    % Near optimality,
    %
    %   fNumerator(x)/(1-x^2) ~= -kappa*u_H'(x)/u_H(x).
    %
    % Hence the integral below reconstructs log(u_H(x)/u_H(0)).
    logGradient = ...
        xi*fNumeratorHandle(xt)./(1-xt.^2);

    logEigenfunction(i) = ...
        -trapz(t,logGradient)/kappa;
end

% The multiplicative normalization is arbitrary. Subtracting the maximum
% before exponentiation prevents overflow without changing the profile.
uSymApprox = exp(logEigenfunction-max(logEigenfunction));

%% Transform the recovered eigenfunction to the plotted density profile

physicalX = L*xPlot;
potential = double_well_potential(physicalX,alpha);

% Multiplication by the Gibbs half-weight converts the recovered
% symmetrized eigenfunction to the quasistationary profile. The nontrivial
% eigenfunction factor comes entirely from the optimized null Lagrangian;
% no exact density is inserted into the computation.
logQsdApprox = logEigenfunction ...
    - 0.5*potential/(sigmaNoise^2/L^2);
qsdApprox = exp(logQsdApprox-max(logQsdApprox));

% Normalize in the physical coordinate so the plotted curve is a density.
qsdApprox = qsdApprox/trapz(physicalX,qsdApprox);

%% Rayleigh quotient upper estimate from the same reconstruction

normU = trapz(xPlot,uSymApprox.^2);

gradFactor = ...
    uSymApprox.*fNumeratorHandle(xPlot)./(1-xPlot.^2);
normGradU = (L^2/sigmaNoise^2)*trapz(xPlot,gradFactor.^2);

wHandle = yalmip_polynomial_to_function(w,x);
normPotential = trapz(xPlot,wHandle(xPlot).*uSymApprox.^2);

lambdaUpper = (normGradU + normPotential)/normU;

%% Print summary

fprintf('\nBounds using degree %d polynomials:\n',degree);
fprintf('  Upper estimate = %.12e\n',lambdaUpper);
fprintf('  Lower bound    = %.12e\n',lambdaLower);
fprintf('  Gap            = %.12e\n',lambdaUpper-lambdaLower);

%% Figure 1: potential landscape

fig1 = figure;
set(fig1,'Color','w','Units','centimeters','Position',[2 2 11 8]);

plot(physicalX,potential, ...
    'Color',S.black, ...
    'LineWidth',3);

xlabel('$\xi$','Interpreter','latex','FontSize',28,'FontWeight','bold');
ylabel('$U(\xi)$','Interpreter','latex','FontSize',28,'FontWeight','bold');

set(gca,'FontSize',22,'TickLabelInterpreter','latex');

box on
grid on
xlim([-L L])

if saveFigures
    export_pdf(fig1,'double_well_potential.pdf');
end

%% Figure 2: reconstructed quasistationary profile

fig2 = figure;
set(fig2,'Color','w');

plot(physicalX,qsdApprox, ...
    'Color',S.blue, ...
    'LineWidth',4);

xlabel('$\xi$','Interpreter','latex','FontSize',36,'FontWeight','bold');
ylabel('$\rho_{0,d}(\xi)$','Interpreter','latex','FontSize',36, ...
    'FontWeight','bold');

set(gca,'FontSize',24,'TickLabelInterpreter','latex');

box on
grid on
xlim([-L L])

if saveFigures
    export_pdf(fig2,'double_well_qsd.pdf');
end

fprintf('\nFinished double-well escape-rate example.\n');

%% Local functions

function U = double_well_potential(x,alpha)

    U = (1/12)*x.^2.*( ...
        -alpha*(6-4*x) ...
        + x.*(3*x-4));
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
