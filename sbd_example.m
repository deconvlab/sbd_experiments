%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('initpkg.m');

%% Generate some synthetic data, activation map values are {0,1}.
% Kernel
p = ceil(10^3);               % Size of the (short) kernel

a0 = randn(p,1);
%a0 = normpdf(1:p, (p+1)/2, p/10)';
a0 = a0/norm(a0);

% Activation / observation
m = 1e2*p;                      % Observation size
theta = 10^-2.5;                % Bernoulli (sparsity) coefficient
dist = @(m,n) randn(m,n);       % Distribution of activations

x0 = (rand(m,1) <= theta) .* dist(m,1);
y = cconv(a0, x0, m);

% Solver properties
params = struct(...
  'solve_lambdas', 0.2*[1/sqrt(p*theta) 1], ...
  'alph', 0, ...
  'iter_ilim', [1 1e3], ...
  'iter_tol', 1e-3);

%% Initialize solver + run some iterations of iPALM
solver = sbd_dq(y, params);
solver.set_ainit(p);

%profile on;
[solver, stats] = solver.solve();
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
