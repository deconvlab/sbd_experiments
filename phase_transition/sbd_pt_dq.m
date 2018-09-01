%% SBD1 Phase transition
clear; %clc; %#ok<*PFBNS>
run('../initpkg.m');

%% Settings
params.theta = 10.^linspace(-2.5, -1.5, 15);
params.p = ceil(10.^linspace(3, 4, 15));
params.trials = 50;
params.saveconvdata = false;

m = @(p) 100*p;
params.xdist = @(theta, p) (rand(m(p),1) <= theta) .* randn(m(p),1);

params.saveconvdata = true;
params.backup = 'results_backup';
params.n_workers = Inf;

lamfac = 0.8;
solparams = @(theta, p) struct('iter_lim', [1 1e3], 'iter_tol', 1e-3, ...
  'solve_lambdas', lamfac *[1/sqrt(p*theta) 1]);
solverfun = @(y, a_init, theta, p) sbd_dq(y, a_init, solparams(theta, p));

%% Run experiment
warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
results = ptloop(solverfun, params);
warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
save('./sbd_pt_dq/tmp1.mat', 'results');
%plotscript;
beep
