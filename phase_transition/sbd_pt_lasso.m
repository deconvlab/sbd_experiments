%% SBD1 Phase transition
clear; %clc; %#ok<*PFBNS>
run('../initpkg.m');

%% Settings
% Data *params
dist = @(m,n) randn(m,n);       % Distribution of activations

thetas = 10.^linspace(-1, -0.6, 10);
p0s = ceil(10.^linspace(1.6, 2, 10));

% Experimental settings
trials = 1;                    % Number of trials
maxit = 1e3;                   % Max iter. & tol. for solver
tol = 1e-3;

%% Containers
tmp = [numel(thetas) numel(p0s)];       % *params
obj = NaN(prod(tmp), trials);
its = NaN(prod(tmp), trials);
times = NaN(trials,1);
idx0 = 0;

%% Run experiment
addpath('./sbd_pt_lasso')
run('../initpkg.m');

loopscript;  save('./sbd_pt_lasso/dq_tmp.mat');
plotscript;  % export_fig ./sbd_pt_lasso/tmp.pdf -transparent

rmpath('./sbd_pt_lasso')
beep
