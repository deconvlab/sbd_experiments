clear; clc;
ipalm_path  = '../../sbd-ipalm_v2/';
[script_path, ~, ~] = fileparts(mfilename('fullpath'));
cd(ipalm_path);   initpkg;
cd(script_path);  initpkg;
addpath('./realdata_ldaupdate');

%% Load data
load('realdata.mat');          % Load data and choose p0
Y = Y/max(abs(Y(:)));
resz = 2;
Y = imresize(Y,resz);
p = [10 10] * resz;

%% Initialize an iPALM iterator for the SBD problem.
xpos = true;                % Recover a nonnegative activation map
getbias = true;             % Recover a constant bias
lambda = 1e-1;              % Sparsity regularization parameter

eta = 1.01;                 % Rate to increase lambda
ldamax = 0.5;               % Max lambda

maxit = 6e2;
trials = 20;

s1 = cell(trials,1);  s2 = s1;
for t = 1:trials
    s1{t} = mkcdl({Y}, p, 1, lambda, true, true);
    s2{t} = copy(s1{t});              % Accelerated
end

updates = [ 1 2:2:10 ...    % when to print updates
            50:50:200 ...
            400:200:max(maxit)];
updates = 1:maxit;

%% Run iterations
A1s = repmat({NaN(p)}, [trials maxit]);  A2s = A1s;
i0 = 1;  loopscript;
save('./realdata_ldaupdate/tmp.mat');

dotscript;
clearvars A1s A2s
save('./realdata_ldaupdate/rd_tmp.mat')

%% Plot dot product
load('./realdata_ldaupdate/rd_tmp.mat')
subplot(212);  plot([mean(dot1,1); mean(dot2,1)]');

%% Finished experiment
rmpath('./realdata_ldaupdate');
beep;
