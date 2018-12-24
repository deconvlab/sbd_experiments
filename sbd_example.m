%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('initpkg.m');

%% Generate some synthetic data.
% Params
p = 500;               % Size of the (short) kernel
m = 1e2*p;                      % Observation size
theta = 10^-3.5;                % Bernoulli (sparsity) coefficient
dist = @(m,n) randn(m,n);       % Distribution of activations

a0 = randn(p,1);
%a0 = normpdf(1:p, (p+1)/2, p/10)';
a0 = a0/norm(a0);
x0 = (rand(m,1) <= theta) .* dist(m,1);
y = cconv(a0, x0, m);

%% Initialize solver + run some iterations of iPALM
% Solver properties
params = struct(...
  'solve_lambdas', 0.2*[1/sqrt(p*theta)], ...
  'alph', 0.9, ...
  'iter_ilim', [1 1e2], ...
  'iter_tol', 1e-3, ...
  'backtrack', NaN, ...
  'refine_iters', [] ...
);

params = struct('alph', 0, 'refine_iters', []);

solver = sbd_dq(y, params);
solver.set_ainit(solver.data_init(p,1));

%profile on;
[solver, stats] = solver.solve();
%profile off; profile viewer
stats %#ok<NOPTS>

% Plot results
subplot(141); plot([y cconv(solver.a, solver.x, m)]); xlim([1 m]);
subplot(142); plot(a0); hold on;
plot(solver.a); hold off; xlim([1 max(numel(a0), numel(solver.a))]);
subplot(143); stem([x0 solver.x], '.'); xlim([1 m]);
subplot(144); plot(cell2mat({stats.costs}));
maxdotshift(a0, solver.a, 0)

%% Done
beep
