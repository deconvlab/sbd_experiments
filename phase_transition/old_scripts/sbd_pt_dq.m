%% SBD1 Phase transition
clear; %clc; %#ok<*PFBNS>


%% Settings
nticks_thetap = [5 5];
trials = 5;                     % Number of trials
save_kernels = false;

maxit = 1e3;                    % Max iter. & tol. for solver
tol = 1e-3;

%thetas = 10.^linspace(-2.5, -1.5, 5);
%p0s = ceil(10.^linspace(2.5, 4.5, 5));
thetas = 10.^linspace(-2.5, -2, nticks_thetap(1));
p0s = ceil(10.^linspace(3, 4, nticks_thetap(2)));

dist = @(m,n) randn(m,n);       % Distribution of activations


%% Containers
tmp = [numel(thetas) numel(p0s)];       % *params
obj = NaN(prod(tmp), trials);
its = NaN(prod(tmp), trials);
times = NaN(trials,1);
idx0 = 0;

%% Run experiment
addpath('./sbd_pt_dq')
run('../initpkg.m');

solverfun = @sbd_dq;
loopscript;  save('./sbd_pt_dq/tmp.mat');
plotscript;

rmpath('./sbd_pt_dq')
beep
