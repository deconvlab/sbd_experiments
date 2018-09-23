%% SBD1 Phase transition
clear; %clc; %#ok<*PFBNS>
run('../initpkg.m');

%% Settings
params.vars = { ...
  'log_t'   linspace(-2.5, -2.0, 10); ...
  'log_p'   linspace(3.2, 4.2, 10) ...
};
params.trials = 10;

params.saveconvdata = false;
params.backup = 'results_backup';
params.n_workers = Inf;

m = @(log_p) 100*round(10^log_p);
params.xdist = @(log_t, log_p) (rand(m(log_p),1) <= 10^log_t) .* randn(m(log_p),1);
params.adist = @(~, log_p) randn(round(10^log_p), 1);

lamfac = 0.8;
solparams = @(theta, p) struct('iter_lim', [1 1e3], 'iter_tol', 1e-3, ...
  'solve_lambdas', lamfac *[1/sqrt(p*theta) 1]);
solverfun = @(y, a_init, log_t, log_p) sbd_dq(y, a_init, solparams(10^log_t, round(10^log_p)));

%% Run experiment
warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
results = loop2var(solverfun, params);
warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

results.vars{1,1} = '\mathrm{Sparsity\ rate\ } \log_{10}(\theta)';
results.vars{2,1} = '\mathrm{Signal\ length\ } \log_{10}(p_0)';
save('./sbd_pt_dq/tmp1.mat', 'results');

%plotscript;
beep
