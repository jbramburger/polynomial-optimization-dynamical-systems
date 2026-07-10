%% rational_sos_recovery.m
% Polynomial Optimization Methods for Dynamical Systems
% Auxiliary Functions, Sum-of-Squares Programming, and Applications
%
% Chapter: Extensions and Frontiers
% Example: Exact rational recovery of a numerical SOS certificate
%
% This script illustrates how a floating-point SOS Gram matrix can be
% converted into an exact rational certificate for
%
%     p(x) = x^4 + 2*x^2 + 1.
%
% Using the monomial basis
%
%     z(x) = [1; x; x^2],
%
% every exact Gram matrix satisfying p = z'Qz has the form
%
%             [1      0       t  ]
%     Q(t) =  [0    2-2t      0  ].
%             [t      0       1  ]
%
% The script:
%
%   1. solves a numerical semidefinite feasibility problem for a Gram matrix,
%   2. rounds its independent entries to nearby rational numbers,
%   3. projects the rounded matrix onto the exact coefficient constraints,
%   4. verifies the polynomial identity using exact integer arithmetic, and
%   5. verifies positive semidefiniteness using exact LDL' sign conditions.
%
% Pressing Run reproduces all numerical and exact verification output.
%
% Requirements:
%   - MATLAB
%   - YALMIP
%   - MOSEK
%
% No Symbolic Math Toolbox is required.
%
% -------------------------------------------------------------------------

clear; clc; close all;
yalmip('clear');
format long g;

%% User options

rationalTolerance = 1e-2;

solverName = 'mosek';
verboseSolver = 0;

%% Solver options

opts = sdpsettings( ...
    'solver',solverName, ...
    'verbose',verboseSolver);

%% Numerical Gram-matrix SDP

Q = sdpvar(3,3,'symmetric');

% Coefficient matching for
%
%     [1 x x^2] Q [1; x; x^2]
%         = 1 + 2x^2 + x^4.

constraints = [
    Q(1,1) == 1, ...
    2*Q(1,2) == 0, ...
    2*Q(1,3) + Q(2,2) == 2, ...
    2*Q(2,3) == 0, ...
    Q(3,3) == 1, ...
    Q >= 0
    ];

fprintf('\n============================================================\n');
fprintf('Exact rational recovery of an SOS certificate\n');
fprintf('============================================================\n\n');

fprintf('Polynomial:\n');
fprintf('  p(x) = x^4 + 2x^2 + 1\n\n');

fprintf('Solving numerical Gram-matrix feasibility problem...\n');

diagnostics = optimize(constraints,[],opts);

if diagnostics.problem ~= 0
    error('The SDP solver failed: %s',diagnostics.info);
end

Qnumerical = value(Q);
Qnumerical = 0.5*(Qnumerical + Qnumerical.');

tNumerical = Qnumerical(1,3);
minimumEigenvalue = min(eig(Qnumerical));

fprintf('\nNumerical Gram matrix Qtilde:\n');
disp(Qnumerical);

fprintf('Numerical free parameter ttilde = %.16g\n',tNumerical);
fprintf('Minimum numerical eigenvalue   = %.6e\n',minimumEigenvalue);

%% Entrywise rational approximation

% Rationalize the independent entries of the numerical Gram matrix. This
% entrywise approximation need not satisfy the coefficient identities
% exactly.

[n11,d11] = rat(Qnumerical(1,1),rationalTolerance);
[n13,d13] = rat(Qnumerical(1,3),rationalTolerance);
[n22,d22] = rat(Qnumerical(2,2),rationalTolerance);
[n33,d33] = rat(Qnumerical(3,3),rationalTolerance);

fprintf('\n------------------------------------------------------------\n');
fprintf('Entrywise rational approximation\n');
fprintf('------------------------------------------------------------\n\n');

fprintf('r11 = %d/%d\n',n11,d11);
fprintf('r13 = %d/%d\n',n13,d13);
fprintf('r22 = %d/%d\n',n22,d22);
fprintf('r33 = %d/%d\n',n33,d33);

% This matrix is constructed in floating point only for display.

Qrounded = [
    n11/d11, 0,       n13/d13;
    0,       n22/d22, 0;
    n13/d13, 0,       n33/d33
    ];

fprintf('\nEntrywise rounded matrix:\n');
disp(Qrounded);

% The coefficient of x^2 in z'Qz is 2Q13 + Q22. Independent entrywise
% rationalization generally causes this value to differ slightly from 2.

roundedCoefficientResidual = ...
    2*(n13/d13) + n22/d22 - 2;

fprintf('Rounded x^2 coefficient residual = %.6e\n', ...
    roundedCoefficientResidual);

%% Exact projection onto the affine coefficient constraints

% Project the rounded entries r13 and r22 onto the affine family
%
%     Q(t) = [1, 0, t; 0, 2-2t, 0; t, 0, 1].
%
% Minimizing the Frobenius distance gives
%
%     t = (2 + r13 - r22)/3.
%
% The calculation below is performed using integer fraction arithmetic:
%
%     r13 = n13/d13,
%     r22 = n22/d22.

projectedNumerator = ...
    2*d13*d22 + n13*d22 - n22*d13;

projectedDenominator = 3*d13*d22;

[projectedNumerator,projectedDenominator] = ...
    reduce_fraction(projectedNumerator,projectedDenominator);

n = projectedNumerator;
d = projectedDenominator;

fprintf('\n------------------------------------------------------------\n');
fprintf('Projection onto the exact coefficient constraints\n');
fprintf('------------------------------------------------------------\n\n');

fprintf('Projected rational parameter:\n');
fprintf('  t = %d/%d = %.16g\n',n,d,n/d);

% Floating-point display of the exact rational matrix.

Qrational = [
    1,         0, n/d;
    0, 2-2*n/d,   0;
    n/d,       0,   1
    ];

fprintf('\nProjected rational Gram matrix Qrat:\n');
disp(Qrational);

fprintf('Distance ||Qrat-Qtilde||_F = %.6e\n', ...
    norm(Qrational-Qnumerical,'fro'));

%% Exact polynomial identity verification

% Multiplying Qrational by d produces the integer matrix
%
%             [d       0       n  ]
%     dQrat = [0    2d-2n      0  ].
%             [n       0       d  ]
%
% The coefficients of
%
%     [1 x x^2](dQrat)[1; x; x^2]
%
% can therefore be compared exactly with those of
%
%     d(1 + 2x^2 + x^4).

Qscaled = [
    d,       0,       n;
    0,       2*d-2*n, 0;
    n,       0,       d
    ];

recoveredScaledCoefficients = [
    Qscaled(1,1);
    2*Qscaled(1,2);
    2*Qscaled(1,3) + Qscaled(2,2);
    2*Qscaled(2,3);
    Qscaled(3,3)
    ];

exactScaledCoefficients = d*[1; 0; 2; 0; 1];

identityResidual = ...
    recoveredScaledCoefficients - exactScaledCoefficients;

fprintf('\n------------------------------------------------------------\n');
fprintf('Exact polynomial identity verification\n');
fprintf('------------------------------------------------------------\n\n');

fprintf('Scaled coefficient residual:\n');
disp(identityResidual);

if any(identityResidual ~= 0)
    error('The projected rational matrix does not satisfy the identity exactly.');
end

fprintf('Polynomial identity verified exactly.\n');

%% Exact positive-semidefiniteness verification

% For
%
%             [1      0       t  ]
%     Q(t) =  [0    2-2t      0  ],
%             [t      0       1  ]
%
% the exact LDL' pivots may be written as
%
%     D11 = 1,
%     D22 = 2 - 2t,
%     D33 = 1 - t^2.
%
% With t = n/d and d > 0, their signs are determined by
%
%     1,
%     2d - 2n,
%     d^2 - n^2.

pivot1Numerator = 1;
pivot2Numerator = 2*d - 2*n;
pivot3Numerator = d^2 - n^2;

fprintf('\n------------------------------------------------------------\n');
fprintf('Exact positive-semidefiniteness verification\n');
fprintf('------------------------------------------------------------\n\n');

fprintf('Exact LDL'' sign data:\n');
fprintf('  D11             = %d\n',pivot1Numerator);
fprintf('  d D22           = %d\n',pivot2Numerator);
fprintf('  d^2 D33         = %d\n',pivot3Numerator);

if pivot1Numerator < 0 || ...
        pivot2Numerator < 0 || ...
        pivot3Numerator < 0

    error('The projected rational Gram matrix is not positive semidefinite.');
end

if pivot1Numerator > 0 && ...
        pivot2Numerator > 0 && ...
        pivot3Numerator > 0

    fprintf('\nThe projected rational Gram matrix is positive definite.\n');

else

    fprintf('\nThe projected rational Gram matrix is positive semidefinite.\n');
end

%% Exact certificate summary

fprintf('\n============================================================\n');
fprintf('Recovered exact SOS certificate\n');
fprintf('============================================================\n\n');

fprintf('The recovered parameter is\n\n');
fprintf('  t = %d/%d.\n\n',n,d);

fprintf('The exact Gram matrix is\n\n');
fprintf('             [1       0       %d/%d]\n',n,d);
fprintf('  Qrat   =   [0       %d/%d   0   ]\n',2*d-2*n,d);
fprintf('             [%d/%d   0       1   ].\n\n',n,d);

fprintf('It satisfies the exact identity\n\n');
fprintf('  x^4 + 2x^2 + 1 = [1 x x^2] Qrat [1; x; x^2]\n\n');
fprintf('and Qrat is positive semidefinite.\n');

fprintf('\nFinished rational SOS recovery example.\n');

%% Local function

function [n,d] = reduce_fraction(n,d)
% Reduce the integer fraction n/d and ensure that its denominator is
% positive.

    if d == 0
        error('Cannot reduce a fraction with zero denominator.');
    end

    commonFactor = gcd(abs(n),abs(d));

    n = n/commonFactor;
    d = d/commonFactor;

    if d < 0
        n = -n;
        d = -d;
    end
end