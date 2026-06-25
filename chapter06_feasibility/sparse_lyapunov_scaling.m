%% sparse_lyapunov_scaling_demo.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Sum-of-Squares Certificates for Stability, Invariance, and Control
% Example: Sparse SOS Lyapunov searches
%
% This script demonstrates how sparsity can reduce the size of SOS
% Lyapunov searches for nearest-neighbour polynomial systems.
%
% The script:
%
%   1. Plots the structured Gram matrix sparsity pattern associated with
%      clique-based SOS constraints on a cycle graph.
%   2. Compares dense quartic Lyapunov searches with sparse structured
%      quartic Lyapunov searches for the cyclic system
%
%          x_i' = -x_i - x_i^3
%                 - gamma x_i (x_{i-1}^2 + x_{i+1}^2),
%
%      with cyclic indexing.
%   3. Reports coefficient counts, monomial counts, Gram basis sizes,
%      solver statuses, and solve times.
%
% Pressing Run should reproduce the numerical output and save the figure
%
%     structured_gram_sparsity.pdf
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

solver_name      = 'mosek';
gamma            = 0.15;
eps_pos          = 1e-4;
eps_decay        = 1e-4;

n_graph          = 12;
n_dense_max      = 12;
n_sparse_max     = 30;
dense_basis_cap  = 600;

verbose_solver   = 0;
makePlots        = true;
saveFigures      = true;

opts = sdpsettings('solver',solver_name,'verbose',verbose_solver);

%% Book plotting style

purple = [170, 90, 160]/255;

%% Part I: structured Gram matrix sparsity pattern

fprintf('\n============================================================\n');
fprintf('Structured Gram matrix sparsity pattern\n');
fprintf('============================================================\n\n');

fprintf('Cycle graph with n = %d vertices.\n', n_graph);
fprintf('Each clique contains variables (x_i,x_{i+1}).\n');

% Degree-6 SOS constraint in two variables uses monomials of degree <= 3.
Nlocal  = nchoosek(2+3,3);
Nstruct = n_graph*Nlocal;

fprintf('Local monomial basis size: %d\n', Nlocal);
fprintf('Structured Gram matrix size: %d x %d\n\n', Nstruct, Nstruct);

S = sparse(Nstruct,Nstruct);

for i = 1:n_graph
    idx = (i-1)*Nlocal + (1:Nlocal);
    S(idx,idx) = 1;
end

fig0 = figure;
set(fig0,'Units','centimeters','Position',[2 2 10 10]);

spy(S,6)
h = findobj(gca,'Type','Line');
set(h,'Color',purple)

axis ij
axis square
axis tight
grid on
box on

set(gca,'FontSize',24,'TickLength',[0 0])
xlabel('row','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('column','Interpreter','latex','FontSize',24,'FontWeight','bold');

if saveFigures
    save_pdf_figure(fig0,'structured_gram_sparsity.pdf');
end

%% Part II: dense quartic SOS searches

dense_results  = [];
sparse_results = [];

fprintf('\n============================================================\n');
fprintf('Sparse Lyapunov scaling demo\n');
fprintf('============================================================\n\n');

fprintf('gamma     = %.4f\n', gamma);
fprintf('eps_pos   = %.1e\n', eps_pos);
fprintf('eps_decay = %.1e\n\n', eps_decay);

fprintf('------------------------------------------------------------\n');
fprintf('Dense full quartic searches\n');
fprintf('------------------------------------------------------------\n');

for n = 2:n_dense_max
    fprintf('\nDense search: n = %d\n',n);

    yalmip('clear');
    x = sdpvar(n,1);

    xp = x([2:n,1]);
    xm = x([n,1:n-1]);

    f = -x - x.^3 - gamma*x.*(xm.^2 + xp.^2);

    [Vd,cd] = polynomial(x,4);
    Vddot = jacobian(Vd,x)*f;

    p1 = sum(x.^2);
    p2 = sum(x.^2);

    [~,monV]    = coefficients(Vd,x);
    [~,monVdot] = coefficients(clean(-Vddot,1e-12),x);

    ncoeff_dense = length(cd);
    nmon_V       = length(monV);
    nmon_Vdot    = length(monVdot);

    gram_deg4 = nchoosek(n+2,2);
    gram_deg6 = nchoosek(n+3,3);

    fprintf('  coefficients in dense V        : %d\n', ncoeff_dense);
    fprintf('  monomials in V                 : %d\n', nmon_V);
    fprintf('  monomials in -Vdot             : %d\n', nmon_Vdot);
    fprintf('  Gram basis size, degree 4 SOS  : %d\n', gram_deg4);
    fprintf('  Gram basis size, degree 6 SOS  : %d\n', gram_deg6);

    if gram_deg6 > dense_basis_cap
        fprintf('  STOP: dense degree-6 Gram basis exceeds cap (%d).\n', ...
            dense_basis_cap);
        break
    end

    Fd = [
        sos(Vd - eps_pos*p1), ...
        sos(-Vddot - eps_decay*p2), ...
        replace(Vd,x,zeros(n,1)) == 0, ...
        replace(Vd,x,[1; zeros(n-1,1)]) == 1
        ];

    tic
    sold = solvesos(Fd,[],opts,cd);
    td = toc;

    fprintf('  solver status code             : %d\n', sold.problem);
    fprintf('  solver message                 : %s\n', sold.info);
    fprintf('  elapsed time (s)               : %.3f\n', td);

    dense_results = [dense_results; ...
        n, ncoeff_dense, nmon_V, nmon_Vdot, gram_deg4, gram_deg6, ...
        td, sold.problem]; 

    if sold.problem ~= 0
        fprintf('  Dense search failed here. Stopping dense sweep.\n');
        break
    end
end

%% Part III: sparse structured quartic SOS searches

fprintf('\n\n------------------------------------------------------------\n');
fprintf('Sparse structured quartic searches\n');
fprintf('------------------------------------------------------------\n');

for n = 2:n_sparse_max
    fprintf('\nSparse structured search: n = %d\n',n);

    yalmip('clear');
    x = sdpvar(n,1);

    xp = x([2:n,1]);
    xm = x([n,1:n-1]);

    f = -x - x.^3 - gamma*x.*(xm.^2 + xp.^2);

    sdpvar a b c
    Vs = a*sum(x.^2) + b*sum(x.^4) + c*sum(x.^2.*xp.^2);
    Vsdot = jacobian(Vs,x)*f;

    p1 = sum(x.^2);
    p2 = sum(x.^2);

    [~,monV]    = coefficients(Vs,x);
    [~,monVdot] = coefficients(clean(-Vsdot,1e-12),x);

    ncoeff_sparse = 3;
    nmon_V        = length(monV);
    nmon_Vdot     = length(monVdot);

    gram_deg4 = nchoosek(n+2,2);
    gram_deg6 = nchoosek(n+3,3);

    fprintf('  coefficients in sparse V       : %d\n', ncoeff_sparse);
    fprintf('  monomials in V                 : %d\n', nmon_V);
    fprintf('  monomials in -Vdot             : %d\n', nmon_Vdot);
    fprintf('  Gram basis size, degree 4 SOS  : %d\n', gram_deg4);
    fprintf('  Gram basis size, degree 6 SOS  : %d\n', gram_deg6);

    Fs = [
        sos(Vs - eps_pos*p1), ...
        sos(-Vsdot - eps_decay*p2), ...
        replace(Vs,x,zeros(n,1)) == 0, ...
        replace(Vs,x,[1; zeros(n-1,1)]) == 1
        ];

    tic
    sols = solvesos(Fs,[],opts,[a;b;c]);
    ts = toc;

    fprintf('  solver status code             : %d\n', sols.problem);
    fprintf('  solver message                 : %s\n', sols.info);
    fprintf('  elapsed time (s)               : %.3f\n', ts);

    sparse_results = [sparse_results; ...
        n, ncoeff_sparse, nmon_V, nmon_Vdot, gram_deg4, gram_deg6, ...
        ts, sols.problem]; 

    if sols.problem == 0
        fprintf('  recovered coefficients         : a = %.6g, b = %.6g, c = %.6g\n', ...
            value(a), value(b), value(c));
    else
        fprintf('  Sparse structured search failed here. Stopping sparse sweep.\n');
        break
    end
end

%% Part IV: summary

fprintf('\n\n============================================================\n');
fprintf('Summary\n');
fprintf('============================================================\n');

if ~isempty(dense_results)
    fprintf('\nDense results columns:\n');
    fprintf('[n, #coeffs, #mon(V), #mon(-Vdot), Gram4, Gram6, time, status]\n');
    disp(dense_results);

    dense_success = dense_results(dense_results(:,8)==0,1);
    if isempty(dense_success)
        largest_dense = NaN;
    else
        largest_dense = max(dense_success);
    end

    fprintf('Largest n solved by dense full quartic search: %g\n', ...
        largest_dense);
end

if ~isempty(sparse_results)
    fprintf('\nSparse results columns:\n');
    fprintf('[n, #coeffs, #mon(V), #mon(-Vdot), Gram4, Gram6, time, status]\n');
    disp(sparse_results);

    sparse_success = sparse_results(sparse_results(:,8)==0,1);
    if isempty(sparse_success)
        largest_sparse = NaN;
    else
        largest_sparse = max(sparse_success);
    end

    fprintf('Largest n solved by sparse structured search: %g\n', ...
        largest_sparse);
end

%% Part V: optional scaling plots

if makePlots
    fig1 = figure;
    set(fig1,'Units','centimeters','Position',[2 2 13 9]);
    hold on

    if ~isempty(dense_results)
        idx = dense_results(:,8)==0;
        plot(dense_results(idx,1), dense_results(idx,7), ...
            'o-', 'LineWidth', 2);
    end

    if ~isempty(sparse_results)
        idx = sparse_results(:,8)==0;
        plot(sparse_results(idx,1), sparse_results(idx,7), ...
            's-', 'LineWidth', 2);
    end

    set(gca,'FontSize',18)
    xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','bold')
    ylabel('solve time (s)','Interpreter','latex','FontSize',20,'FontWeight','bold')
    legend({'dense full quartic','sparse structured quartic'}, ...
        'Interpreter','latex','Location','northwest')
    grid on
    box on

    if saveFigures
        save_pdf_figure(fig1,'sparse_lyapunov_scaling_times.pdf');
    end

    fig2 = figure;
    set(fig2,'Units','centimeters','Position',[2 2 13 9]);
    hold on

    if ~isempty(dense_results)
        plot(dense_results(:,1), dense_results(:,2), ...
            'o-', 'LineWidth', 2);
    end

    if ~isempty(sparse_results)
        plot(sparse_results(:,1), sparse_results(:,2), ...
            's-', 'LineWidth', 2);
    end

    set(gca,'FontSize',18)
    xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','bold')
    ylabel('number of decision coefficients', ...
        'Interpreter','latex','FontSize',20,'FontWeight','bold')
    legend({'dense full quartic','sparse structured quartic'}, ...
        'Interpreter','latex','Location','northwest')
    grid on
    box on

    if saveFigures
        save_pdf_figure(fig2,'sparse_lyapunov_scaling_coeffs.pdf');
    end
end

fprintf('\nFinished sparse Lyapunov scaling demo.\n');

%% Local function

function save_pdf_figure(fig,fileName)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n', fileName);
end