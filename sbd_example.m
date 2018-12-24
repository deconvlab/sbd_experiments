%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('initpkg.m');

%% Generate some synthetic data.
% Params
p = 5e2;                      % Size of the (short) kernel
m = 1e3*p;                    % Observation size
theta = 2/p;                  % Bernoulli (sparsity) coefficient
dist = @(m) randn(m,1);       % Distribution of activations

a0 = randn(p,1);
%a0 = normpdf(1:p, (p+1)/2, p/10)';
a0 = a0/norm(a0);
x0 = (rand(m,1) <= theta) .* dist(m);
y = cconv(a0, x0, m);

%% Initialize solver + run some iterations of iPALM
% Solver properties
params = struct();
params.solve_lambdas = 2e-1*[1/sqrt(p*theta) 1];
params.alph = 0;
params.iter_lim = [1 1e3];
params.iter_tol = 1e-3;
params.backtrack = [0.1 0.1];
params.refine_iters = [];

solver = sbd_lasso(params);
solver = sbd_dq(params);
solver.y = y;
solver.a = data_ainit(solver, p,1);

%profile on;
[solver, stats] = solver.solve();
%profile off; profile viewer
clc; stats(1).it

% Plot results
subplot(141); plot([y cconv(solver.a, solver.x, m)]); xlim([1 m]);
subplot(142); plot(a0); hold on;
plot(solver.a); hold off; xlim([1 max(numel(a0), numel(solver.a))]);
subplot(143); stem([x0 solver.x], '.'); xlim([1 m]);
subplot(144); plot(cell2mat({stats.costs}));
maxdotshift(a0, solver.a, 0)

%% Done
beep
