%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('initpkg.m');

%% Generate some synthetic data, activation map values are {0,1}.
% Kernel
p = 20;                         % Size of the (short) kernel

a0 = randn(p,1);
%a0 = normpdf(1:p, (p+1)/2, p/10)';
a0 = a0/norm(a0);

% Activation / observation
m = 1024*5;                     % Observation size
theta = 1e-1;                   % Bernoulli (sparsity) coefficient
dist = @(m,n) randn(m,n);       % Distribution of activations

x0 = (rand(m,1) <= theta) .* dist(m,1);
y = cconv(a0, x0, m);

% Solver properties
lambda = 1e-1;                  % Sparsity regularization parameter
maxit = 2e2;                    % iterations in initial iPALM solve
tol = 1e-4;

%% Initialize solver + run some iterations of iPALM
solver = sbd_lasso(y, 3*p-2, struct('weights', lambda));
%[solver, stats] = iterate(solver, [10 maxit], tol);
[solver, stats] = solve(solver, [10 maxit], tol, 3);

%% Plot results
subplot(221); plot([y cconv(solver.a, solver.x, m)]); xlim([1 m]);
subplot(222); plot(a0); hold on;
plot(solver.a); hold off; xlim([1 max(numel(a0), numel(solver.a))]);
subplot(223); stem([x0 solver.x], '.'); xlim([1 m]);
subplot(224); plot(cell2mat({stats.costs}));
maxdotshift(a0, solver.a, 0)

