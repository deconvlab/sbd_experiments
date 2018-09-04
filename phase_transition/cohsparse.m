%% SBD1 Phase transition
clear; %clc; %#ok<*PFBNS>
run('../initpkg.m');

%% Settings
params.trials = 5;
params.saveconvdata = false;
params.backup = 'results_backup';
params.n_workers = 0; %Inf;

log_p = 3;
p = round(10^log_p);
m = p * 100;

params.vars = { ...
  'log_t'     linspace(-2.5, -2, 5); ...
  'log_pmu'   linspace(log_p-1, log_p, 5) ...
};

params.xdist = @(log_t, ~) (rand(m,1) <= 10^log_t) .* randn(m,1);
params.adist = @(~, log_pmu) [randn(round(10^log_pmu),1); zeros(p-round(10^log_pmu),1)];

lamfac = 0.8;
solparams = @(theta, p) struct('iter_lim', [1 1e3], 'iter_tol', 1e-3, ...
  'solve_lambdas', lamfac *[1/sqrt(p*theta) 1]);
solverfun = @(y, a_init, log_t, ~) sbd_dq(y, a_init, solparams(10^log_t, p));

%% Run experiment
warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
results = loop2var(solverfun, params);
warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

%% Save results
results.vars{1,1} = '\log_{10}(\theta)';
results.vars{2,1} = '\log_{10}(p_\mu)';
save('./cohsparse/tmp.mat', 'results');

%plotscript;
beep
