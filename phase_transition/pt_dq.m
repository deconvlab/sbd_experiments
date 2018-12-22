%% SBD1 Phase transition
clear; %clc; %#ok<*PFBNS>
run('../initpkg.m');

%% Settings
ngrid = 3;              % Number of elements per grid
params.trials = 2;      % Number of trials per choice of parameters

params.vars = { ...
  %'log_t'   linspace(-3.5, -3, ngrid); ...   % log theta
  'log_t'   linspace(-4, -3, ngrid); ...

  %'log_p'   linspace(3.5, 4.5, ngrid); ...   % log p
  'log_p'   linspace(3, 4, ngrid); ...
};

params.saveconvdata = true;        % whether or not to save variables
params.backup = 'results_backup';   % variable to update as the trials go
params.n_workers = Inf;

% How to pick x0, a0, and a_init
m = @(log_p) 100*round(10^log_p);
params.gen.x0 = @(log_t, log_p) (rand(m(log_p),1) <= 10^log_t) .* randn(m(log_p),1);
params.gen.a0 = @(~, log_p) randn(round(10^log_p), 1);
params.gen.ainit = @(~,log_p) randn(round(10^log_p), 1);

% Solver parameters for differing theta and p
% e.g. for 5 refinement steps; 10 fista and 20 agd iterations each
%   refine_iters = ones(5,1)*[10 20];
solparams = @(theta, p) struct( ...
  'iter_lim', [1 1e3], ...      % min / max iterations
  'iter_tol', 1e-3, ...         % tolerence until exiting iterations
  'solve_lambdas', 0.2*[1/sqrt(p*theta) 1], ...   % sequence of lambdas
  'refine_iters', [] ...        % xsolve and agd iterations
);

% How the solver gets initialized as parameters change
solverfun = @(y, a_init, log_t, log_p) ...
  sbd_dq(y, a_init, solparams(10^log_t, round(10^log_p)));

%% Run experiment
warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
results = loop2var(solverfun, params);
warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

% Labels for plots
results.vars{1,1} = '\mathrm{Sparsity\ rate\ } \log_{10}(\theta)';
results.vars{2,1} = '\mathrm{Signal\ length\ } \log_{10}(p_0)';
save('./pt_dq/tmp1.mat', 'results');

%plotscript;
beep
