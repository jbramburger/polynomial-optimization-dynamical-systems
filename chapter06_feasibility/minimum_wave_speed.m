%% minimum_wave_speed.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: SOS bounds on a minimum traveling-wave speed
%
% This script computes SOS upper and lower bounds for the minimum wave
% speed in a non-KPP isothermal diffusion model.
%
% The script combines:
%
%   1. an upper-bound certificate based on a trapping-region formulation,
%   2. a lower-bound certificate based on a volume/auxiliary-function method.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     nonkpp_minimum_wave_speed_bounds.pdf
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

n_reaction = 2;

D_values = 0.02:0.02:2.00;

upper_degree = 12;
lower_degree = 2;

lower_lambda = 1e3;

bisect_tol = 1e-4;

solver_name = 'mosek';
verbose_solver = 0;

saveFigure = true;
figureFile = 'nonkpp_minimum_wave_speed_bounds.pdf';

opts = sdpsettings( ...
    'solver',solver_name, ...
    'verbose',verbose_solver, ...
    'cachesolvers',1);

%% Book plotting style

blue  = [0, 92, 175]/255;
red   = [200, 50, 50]/255;
black = [0, 0, 0];

%% Storage

upper_bound = nan(size(D_values));
lower_bound = nan(size(D_values));

fprintf('\n============================================================\n');
fprintf('SOS bounds on a minimum traveling-wave speed\n');
fprintf('============================================================\n\n');

fprintf('reaction exponent n = %d\n', n_reaction);
fprintf('upper degree        = %d\n', upper_degree);
fprintf('lower degree        = %d\n', lower_degree);
fprintf('lower lambda        = %.3g\n\n', lower_lambda);

%% Main loop over D

for k = 1:length(D_values)

    D = D_values(k);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('D = %.4f\n', D);
    fprintf('------------------------------------------------------------\n');

    %% Upper bound

    if abs(D - 1) < 1e-12

        upper_bound(k) = 1/sqrt(2);
        fprintf('Upper bound: D = 1 endpoint, using c = %.8f\n', upper_bound(k));

    else

        if D < 1
            c_left = 1e-4;
            c_right = 1/sqrt(2);
        else
            c_left = sqrt(D/2);
            c_right = D;
        end

        upper_bound(k) = bisection_upper_bound( ...
            c_left, c_right, bisect_tol, upper_degree, D, n_reaction, opts);

        fprintf('Upper bound: %.8f\n', upper_bound(k));
    end

    %% Lower bound

    if abs(D - 1) < 1e-12

        lower_bound(k) = 1/sqrt(2);
        fprintf('Lower bound: D = 1 endpoint, using c = %.8f\n', lower_bound(k));

    else

        if D < 1
            c_left = 1e-4;
            c_right = 1/sqrt(2);
        else
            c_left = sqrt(D/2);
            c_right = D;
        end

        lower_bound(k) = bisection_lower_bound( ...
            c_left, c_right, bisect_tol, lower_degree, D, ...
            n_reaction, lower_lambda, opts);

        fprintf('Lower bound: %.8f\n', lower_bound(k));
    end
end

%% Print summary

fprintf('\n============================================================\n');
fprintf('Computed speed bounds\n');
fprintf('============================================================\n');

fprintf('\nColumns: D, lower bound, upper bound\n');
disp([D_values(:), lower_bound(:), upper_bound(:)]);

%% Plot

fig = figure;
set(fig,'Color','w','Units','centimeters','Position',[2 2 13 9]);

hold on
box on

plot(D_values, lower_bound, '-', ...
    'Color', blue, ...
    'LineWidth', 3);

plot(D_values, upper_bound, '--', ...
    'Color', red, ...
    'LineWidth', 3);

plot(D_values, sqrt(D_values/2), ':', ...
    'Color', black, ...
    'LineWidth', 2);

xlabel('$D$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('$c$','Interpreter','latex','FontSize',24,'FontWeight','bold');

legend({'lower bound','upper bound','$\sqrt{D/2}$'}, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'FontSize',18, ...
    'Box','on');

set(gca,'FontSize',24)
grid on

if saveFigure
    save_pdf_figure(fig,figureFile);
end

fprintf('\nFinished non-KPP minimum-wave-speed example.\n');

%% Local functions

function c_upper = bisection_upper_bound(c_left, c_right, tol, d, D, n, opts)

    while abs(c_right - c_left) >= tol

        c_mid = 0.5*(c_left + c_right);

        if D < 1
            feasible = upper_certificate_D_less_than_one(c_mid,d,D,n,opts);
        else
            feasible = upper_certificate_D_greater_than_one(c_mid,d,D,n,opts);
        end

        if feasible
            c_right = c_mid;
        else
            c_left = c_mid;
        end
    end

    c_upper = c_right;
end

function c_lower = bisection_lower_bound(c_left, c_right, tol, d, D, n, lambda, opts)

    while abs(c_right - c_left) >= tol

        c_mid = 0.5*(c_left + c_right);

        feasible = lower_certificate(c_mid,d,D,n,lambda,opts);

        if feasible
            c_left = c_mid;
        else
            c_right = c_mid;
        end
    end

    c_lower = c_left;
end

function feasible = upper_certificate_D_less_than_one(c,d,D,n,opts)

    r = 1/D - 1;

    sdpvar x y
    p = sdpvar(2,1);

    [tildeF, cF] = polynomial([x y], d);
    F = (1 - y)*tildeF;

    Fx = jacobian(F,x);

    P = -int( ...
        -F ...
        + ((D/c)^2)*(y + r*x)*(1-y)^n ...
        - Fx*((1 - D)*y - D*x*r), ...
        y, y, 1);

    dmult = d;

    [s11, c11] = polynomial([x y],dmult);
    [s12, c12] = polynomial([x y],dmult);
    [s22, c22] = polynomial([x y],dmult);

    [t11, ct11] = polynomial([x y],dmult);
    [t12, ct12] = polynomial([x y],dmult);
    [t22, ct22] = polynomial([x y],dmult);

    S = [s11, s12; s12, s22];
    T = [t11, t12; t12, t22];

    [s4, c4] = polynomial([x y], dmult);
    [s5, c5] = polynomial([x y], dmult);

    Q = [P, F; F, 2];

    cons = [
        sos(F - x*(1-x)*s4 - y*(1-y)*s5), ...
        sos(p'*(Q - x*(1-x)*S - y*(1-y)*T)*p), ...
        sos(p'*S*p), ...
        sos(p'*T*p), ...
        sos(s4), ...
        sos(s5)
        ];

    decision_vars = [cF; c11; c12; c22; ct11; ct12; ct22; c4; c5];

    sol = solvesos(cons, [], opts, decision_vars);

    feasible = (sol.problem == 0);
end

function feasible = upper_certificate_D_greater_than_one(c,d,D,n,opts)

    r = 1 - 1/D;

    sdpvar x y
    p = sdpvar(2,1);

    [tildeF, cF] = polynomial([x y], d);
    F = (1 - y)*tildeF;

    Fx = jacobian(F,x);

    P = -int( ...
        -F ...
        + ((D/c)^2)*(y - r*x)*(1-y)^n ...
        - Fx*((D-1)*y - D*x*r), ...
        y, y, 1);

    dmult = d;

    [s11, c11] = polynomial([x y],dmult);
    [s12, c12] = polynomial([x y],dmult);
    [s22, c22] = polynomial([x y],dmult);

    [t11, ct11] = polynomial([x y],dmult);
    [t12, ct12] = polynomial([x y],dmult);
    [t22, ct22] = polynomial([x y],dmult);

    S = [s11, s12; s12, s22];
    T = [t11, t12; t12, t22];

    [s4, c4] = polynomial([x y], dmult);
    [s5, c5] = polynomial([x y], dmult);

    Q = [P, F; F, 2];

    cons = [
        sos(F - x*(1-x)*(y - r*x)*(1-y)*s5), ...
        sos(p'*(Q - x*(1-x)*S - (y - r*x)*(1-y)*T)*p), ...
        sos(p'*S*p), ...
        sos(p'*T*p), ...
        sos(s4), ...
        sos(s5)
        ];

    decision_vars = [cF; c11; c12; c22; ct11; ct12; ct22; c4; c5];

    sol = solvesos(cons, [], opts, decision_vars);

    feasible = (sol.problem == 0);
end

function feasible = lower_certificate(c,d,D,n,lambda,opts)

    eps_margin = 1e-4;

    if D < 1
        Tscale = 1e2;
        r = 1/D - 1;
    else
        Tscale = 1;
        r = 1 - 1/D;
    end

    sdpvar u v w

    [V, cV] = polynomial([u v w], [d d d]);

    [s1, c1] = polynomial([u v w], d);
    [s2, c2] = polynomial([u v w], d);
    [s3, c3] = polynomial([u v w], d);

    [s4, c4] = polynomial([u v], d);
    [s5, c5] = polynomial([u v], d);

    Vu = jacobian(V,u);
    Vv = jacobian(V,v);
    Vw = jacobian(V,w);

    Vw0 = replace(V,w,0);

    if D < 1

        F1 = (1-D)*v - D*u*r*Tscale;
        F2 = w;
        F3 = -w + ((D/c)^2)*(v + r*u*Tscale)*(1-v)^n;

        cons = [
            replace(V, [u v w], [0 0 0]) == 0, ...
            sos(lambda*(Vu*F1 + Vv*F2 + Vw*F3) ...
                - V ...
                - u*(1-u*Tscale)*s1 ...
                - v*(1-v)*s2 ...
                - w*s3), ...
            sos(-Vw0 - eps_margin*(u+v) ...
                - u*(1-u*Tscale)*s4 ...
                - v*(1-v)*s5), ...
            sos(s1), ...
            sos(s2), ...
            sos(s3), ...
            sos(s4), ...
            sos(s5)
            ];

    else

        F1 = (D-1)*v - D*r*u;
        F2 = w;
        F3 = -w + ((D/c)^2)*(v-r*u)*(1-v)^n;

        cons = [
            replace(V, [u v w], [0 0 0]) == 0, ...
            sos(lambda*(Vu*F1 + Vv*F2 + Vw*F3) ...
                - V ...
                - u*(1-u)*s1 ...
                - (v-r*u)*(1-v)*s2 ...
                - w*s3), ...
            sos(-Vw0 - eps_margin*(u+v) ...
                - u*(1-u)*s4 ...
                - (v-r*u)*(1-v)*s5), ...
            sos(s1), ...
            sos(s2), ...
            sos(s3), ...
            sos(s4), ...
            sos(s5)
            ];
    end

    decision_vars = [cV; c1; c2; c3; c4; c5];

    sol = solvesos(cons, [], opts, decision_vars);

    feasible = (sol.problem == 0);
end

function save_pdf_figure(fig,fileName)

    set(fig, 'PaperUnits','centimeters');
    set(fig, 'Units','centimeters');

    pos = get(fig,'Position');

    set(fig, 'PaperSize', [pos(3) pos(4)]);
    set(fig, 'PaperPositionMode', 'manual');
    set(fig, 'PaperPosition', [0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n', fileName);
end