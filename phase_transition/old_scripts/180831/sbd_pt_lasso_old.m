%% SBD1 Phase transition
clear; clc; %#ok<*PFBNS>
run('../initpkg.m');

%% Experimental parameters
% Data properties
dist = @(m,n) randn(m,n);       % Distribution of activations
p0 = 1e2;                       % Fixed kernel size of a0
p = 2*p0-2;                   % Size of recovery window
%p = ceil(p0*1.1);

trials = 20;                    % Number of trials
tfacs = linspace(0.01,0.8,15);   % How to space theta (see below)
mfacs = linspace(20,200,15);    % How to space m (see below)
%mfacs = 200;

thetas = tfacs;
ms = ceil(p0 * mfacs);

% Solver properties
lambda = [1 1]*1e-1;            % Manually set lambda for each cont. phase
maxit = 1e3;                    % Max iter. & tol. for solver
tol = 1e-3;

%% Initialize solver + run some iterations of iPALM
tmp = [numel(thetas) numel(ms)];
obj = NaN(prod(tmp), trials);
its = NaN(prod(tmp), trials);
a0 = randn(p0,1);  a0 = a0/norm(a0);        % a0 is FIXED

clc;
warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

for idx = 1:prod(tmp)
    fprintf('Testing %d of %d...\n', idx, prod(tmp));
    [i, j] = ind2sub(tmp, idx);
    theta = thetas(i);  m = ms(j);

    start = tic;
% WHAT HAPPENS IN EACH TRIAL:
parfor trial = 1:trials
    % A) Generate x & y: supp(x) must be >= 1
    xgood = false;
    while ~xgood
        x0 = (rand(m,1) <= theta) .* dist(m,1);
        xgood = sum(x0~=0) >= 1;
    end
    y = cconv(a0, x0, m);
    
    % B) Create solver and run continuation sequence
    solver = sbd_lasso(y, p, struct('lambda', lambda(1))); 
    solver = solve(solver, [10 maxit], tol, lambda);
    
    % C) Record statistics
    obj(idx, trial) = maxdotshift(a0, solver.a);
    its(idx, trial) = solver.it;
end
    fprintf('\b\b\b\b: theta = %.3f, m = %d, mean obj. = %.2f.', ...
        theta, m, mean(obj(idx,:)));
    fprintf(' Time elapsed: %.1fs.\n', toc(start));
end
obj = reshape(obj, [tmp trials]);
its = reshape(its, [tmp trials]);
warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
disp('Done.');

% Plots
tmp = {0.9 'flat'};
i = tfacs;  j = mfacs;

clf;
if min(size(obj,1), size(obj,2)) == 1
    % 1D transition plots
    if size(obj,1) == 1;  tmp{3} = j;  else;  tmp{3} = i;  end
    yyaxis left;   
    plot(tmp{3}, median(squeeze(obj),2));  hold on;
    plot(tmp{3}, mean(squeeze(obj),2));
    plot(tmp{3}, min(squeeze(obj), [], 2));
    hold off;
    ylabel('\rho = max_i <S_i[a0], a>');

    yyaxis right;    plot(tmp{3}, mean(squeeze(obj) >= tmp{1},2));
    ylabel(['P[\rho > ' sprintf('%.2f]', tmp{1})]); 
    xlabel('p\theta');
    
else
    % 2D transition plots
    colormap gray; 
    subplot(131); surf(i, j, mean(obj,3)'); 
    view(2); shading(tmp{2});  title('Mean of \rho');
    
    subplot(132); surf(i, j, median(obj,3)'); 
    view(2); shading(tmp{2});  title('Median of \rho');
    
    subplot(133); surf(i, j, mean(obj>=tmp{1},3)'); 
    view(2); shading(tmp{2});  title('Success prob.');
end

% End of experiment
beep
