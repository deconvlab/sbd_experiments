%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('initpkg.m');

%% Generate some synthetic data, activation map values are {0,1}.
% Kernel
p = 1e2;                         % Size of the (short) kernel

a0 = randn(p,1);
%a0 = normpdf(1:p, (p+1)/2, p/10)';
a0 = a0/norm(a0);

% Activation / observation
m = 1e3*2;                      % Observation size
%m = p*3e1;
theta = 1e-1;                   % Bernoulli (sparsity) coefficient
%theta = 1e0/p;
dist = @(m,n) randn(m,n);       % Distribution of activations

x0 = (rand(m,1) <= theta) .* dist(m,1);
y = cconv(a0, x0, m);

% Solver properties
lambda = ones(2,1) * 1e-1;      % Sparsity regularization parameter
maxit = 1e3;                    % iterations in initial iPALM solve
tol = 1e-3;

%% Initialize solver + run some iterations of iPALM
solver = sbd_lasso(y, 3*p-2, struct('lambda', lambda(1)));
%solver = sbd_dq(y, 3*p-2, struct('alph', 0.9, 'lambda', lambda(1)));

%profile on;
[solver, stats] = solve(solver, [10 maxit], tol, lambda);
%profile off; profile viewer

%% Plot results
subplot(141); plot([y cconv(solver.a, solver.x, m)]); xlim([1 m]);
subplot(142); plot(a0); hold on;
plot(solver.a); hold off; xlim([1 max(numel(a0), numel(solver.a))]);
subplot(143); stem([x0 solver.x], '.'); xlim([1 m]);
subplot(144); plot(cell2mat({stats.costs}));
maxdotshift(a0, solver.a, 0)

%% Done
beep
