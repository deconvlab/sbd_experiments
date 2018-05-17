%% SBD1 Phase transition
clear; %clc; %#ok<*PFBNS>

%% Settings
% Data *params
dist = @(m,n) randn(m,n);       % Distribution of activations
%thetas = 10.^linspace(-2.5, -1.5, 5);
%p0s = ceil(10.^linspace(2.5, 4.5, 5));

thetas = 10.^linspace(-2.5, -2, 15);
p0s = ceil(10.^linspace(3, 4, 15));

% Experimental settings
trials = 50;                    % Number of trials
maxit = 1e3;                    % Max iter. & tol. for solver
tol = 1e-3;

%% Containers
tmp = [numel(thetas) numel(p0s)];       % *params
obj = NaN(prod(tmp), trials);
its = NaN(prod(tmp), trials);
times = NaN(trials,1);
idx0 = 0;

%% Run experiment
addpath('./sbd_pt_dq')
run('../initpkg.m');

loopscript;
plotscript;

rmpath('./sbd_pt_dq')
beep
