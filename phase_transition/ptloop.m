function results = ptloop(solverfun, params, results)
%PTLOOP  Run phase transition experiment
%
% RESULTS = PTLOOP(SOLVER, PARAMS, RESULTS)
%
% Arguments:
%   SOLVERFUN:  Supply a function handle to a SBD solver as follows:
%     SOLVER = SOLVERFUN(Y, A_INIT, THETA, P).
%
%   PARAMS:  A struct with the following fields:
%     theta, p:  a 1d real array for the sparsity rate and kernel size.
%
%     xdist:  A function handle for the distribution of x, as follows:
%       X = XDIST(THETA, P).
%
%     trials:  The number of trials to solve for each theta-p
%
%     saveconvdata: Optional bool value. Returns A0, X0 and A for each
%       trial sample if true.
%
%     backup:  Optional string. If nonempty, backs up a copy of the
%       results at the end of each trial using the string as filename.
%       See RESULTS in the Returns section.
%
%     exp_init:  Optional experiment index to start from. Default is 1.
%
%   RESULTS:  A copy of the results to pick off from if the loop was
%     terminated prematurely. See RESULTS in the Returns section.
%
% Returns:
%   RESULTS:  A struct containing all the results from the experiment.
%

  % Pass parameters:
  theta = params.theta;
  p = params.p;
  n_exps = numel(theta) * numel(p);
  trials = params.trials;

  saveconvdata = isfield(params, 'saveconvdata') && params.saveconvdata

  if isfield(params, 'backup')
    backup = params.backup;
  else
    backup = '';
  end

  if isfield(params, 'exp_init')
    exp_init = params.exp_init;
  else
    exp_init = 1;
  end

  if nargin < 3 || ~isstruct(results)
    results.n_exps = n_exps;
    results.theta = theta;
    results.p = p;
    results.trials = trials;

    results.obj = NaN(n_exps, trials);
    results.its = NaN(n_exps, trials);
    results.times = NaN(trials,1);

    if saveconvdata
      results.a0 = cell(n_exps, trials);
      results.x0 = cell(n_exps, trials);
      results.a = cell(n_exps, trials);
    end
  end

  %warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
  for idx = exp_init:n_exps
    fprintf('[%5d|%5d]:  ', idx-1, n_exps-1);
    [i, j] = ind2sub([numel(theta) numel(p)], idx);

    thetai = theta(i); pj = p(j);
    a0 = randn(pj,1);  a0 = a0/norm(a0);

    % WHAT HAPPENS IN EACH TRIAL:
    t0 = tic;
    for trial = 1:trials
      % A) Generate x & y: supp(x) must be >= 1
      xgood = false;
      while ~xgood
        x0 = params.xdist(thetai, pj);
        xgood = sum(x0~=0) >= 1;
      end
      y = cconv(a0, x0, m);

      % B) Create solver and run continuation sequence
      solver = solverfun(y, a0, thetai, pj);
      solver.solve();

      % C) Record statistics
      results.obj(idx, trial) = maxdotshift(a0, solver.a);
      results.its(idx, trial) = solver.it;

      if saveconvdata
        results.a0{idx, trial} = a0;
        results.x0{idx, trial} = sparse(x0);
        results.a{idx, trial} = solver.a;
      end
    end
    results.times(idx) = toc(t0);

    % Print & back up results
    fprintf('theta = %.2E, p0 = %d, mean obj. = %.2f.', ...
      pj, thetai, mean(results.obj(idx,:)));
    fprintf(' Time elapsed: %.1fs.\n', results.times(idx));
    if i == numel(thetas);  disp(' ');  end

    if numel(backup);  assignin('base', backup, results);  end
  end

  %% Final adjustments to results.
  %warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
  results.obj = reshape(results.obj, [tmp trials]);
  results.its = reshape(results.its, [tmp trials]);
  disp('Done.');
end

%#ok<*SAGROW>