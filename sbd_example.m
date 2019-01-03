%% Example - how to use the iPALM code for a CDL problem
clear; clc; 
run('initpkg.m');

%% Generate some synthetic data.
% Params
p = 1e3;                      % Size of the (short) kernel
m = 1e5;                    % Observation size
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
params.solve_lambdas = 5e-1*[1/sqrt(p*theta) 1];    % sequence of lambdas
params.alph = 0;                  % momentum
params.iter_lim = [1 1e3];        % min/max iterations per lambda
params.iter_tol = 1e-3;           % stepsize*grad together must be small
params.backtrack = [0.1 0.1];     % [btdec btslack]; set empty [] to turn off

params.refine_iters = ones(5,1)*[10 20];            % see below
% n_refine x 2 array; set empty [] to turn off
% e.g. for 5 refinement steps; 10 fista and 20 agd iterations each
%   refine_iters = ones(5,1)*[10 20];

solver = sbd_lasso(params);
%solver = sbd_dq(params);
solver.y = y;
solver.a = data_ainit(solver, p,1);

%profile on;
[solver, stats] = solver.solve();
%profile off; profile viewer

%% Plot results
lgd = {'Truth', 'Recovered'};

subplot(411); plot([y cconv(solver.a, solver.x, m)], 'LineWidth', 1.2); 
xlim([1 m]); ylim([1 1.5].*ylim);  title('Observation y'); legend(lgd);

subplot(412); plot(a0, 'LineWidth', 1.2);  hold on;
plot(solver.a, 'LineWidth', 1.2);  hold off; 
xlim([1 max(numel(a0), numel(solver.a))]);  ylim([-1 1]*max(abs(ylim)));  legend(lgd);
title(sprintf('Kernel a:    max_i |<s_i[a_0], a>| = %f', maxdotshift(a0, solver.a, 0)));

subplot(413); stem([x0 solver.x], '.', 'LineWidth', 1.2, 'MarkerSize', 10); 
xlim([1 m]); ylim([1 1.5].*ylim); title('Activation x'); legend(lgd);

subplot(414); 
it = 0;
for i = 1:numel(stats)
  plot(it+(1:stats(i).it), stats(i).costs, '.-', 'LineWidth', 1.2, 'MarkerSize', 10); hold on;
  it = it + stats(i).it;
end
hold off;
xlim([1 it]);
xlabel('Iteration #'); title('Objective');

%% Done
beep
