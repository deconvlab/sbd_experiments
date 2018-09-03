function results = loop2var(solverfun, params, results)
%LOOP2VAR  Run experiment on two variables
%
% RESULTS = LOOP2VAR(SOLVERFUN, PARAMS, RESULTS)
%
% Arguments:
%   SOLVERFUN: Function handle.  When called creates a SBD solver as follows:
%     SOLVER = SOLVERFUN(Y, A_INIT, VAR1, VAR2).
%
%   PARAMS: Struct.  With the following fields:
%     vars: 2x2 cell array.  Each row contains the {VARNAME, VARARRAY}
%       elements, where VARNAME is a string and VARARRAY is a 1d real array.
%
%     xdist: Function handle.  For the distribution of x, as follows:
%       X = XDIST(VAR1, VAR2).
%
%     adist: Function handle.  For the distribution of a, as follows:
%       A = ADIST(VAR1, VAR2). The resulting kernel is then normalized.
%
%     trials: Int.  The number of trials to solve for each theta-p
%
%     saveconvdata (optional): Bool value.  Returns A0, X0 and A for each
%       trial sample if true. Default is false.
%
%     backup (optional): String.  If nonempty, backs up a copy of the results
%       at the end of each experiment to base workspace using the string as
%       variable name. See RESULTS in the Returns section.
%
%     exp_init (optional): Int.  Experiment index to start from. Default is 1.
%
%     n_workers (optional): Int.  Number of workers for parfor. Default is 0.
%
%   RESULTS (optional): Struct.  Results to pick off from if the loop was
%     terminated prematurely. See RESULTS in the Returns section.
%
% Returns:
%   RESULTS: Struct.  Contains all the results from the experiment.
%

  % Pass parameters:
  var1 = params.vars{1,2};
  var2 = params.vars{2,2};
  n_exps = numel(var1) * numel(var2);
  trials = params.trials;
  xdist = params.xdist;
  adist = params.adist;

  saveconvdata = isfield(params, 'saveconvdata') && params.saveconvdata;

  optional_params = { ...
    'saveconvdata', 'false'; ...
    'backup', ''; ...
    'exp_init', 1; ...
    'n_workers', 0 ...
  };

  for idx = 1:size(optional_params,1)
    if isfield(params, optional_params{idx,1})
      eval([optional_params{idx,1} '= params.' optional_params{idx,1} ';']);
    else
      eval([optional_params{idx,1} '= optional_params{idx,2};']);
    end
  end

  if nargin < 3 || ~isstruct(results)
    results.n_exps = n_exps;
    results.vars = vars;
    results.trials = trials;

    results.obj = NaN(n_exps, trials);
    results.obj_ = NaN(n_exps, trials);
    results.its = NaN(n_exps, trials);
    results.times = NaN(n_exps, trials);

    if saveconvdata
      results.a0 = cell(n_exps, trials);
      results.x0 = cell(n_exps, trials);
      results.a = cell(n_exps, trials);
    end
  end

  for idx = exp_init:n_exps
    fprintf('[%5d|%5d]:  ', idx-1, n_exps-1);

    [i, j] = ind2sub([numel(theta) numel(p)], idx);
    v1i = var1(i); v2j = var2(j);

    % Create parfor containers:
    % In ptplot.m, the final objective (obj) will be compared to the
    % objective in the last stats container produced by the solve method
    % of the sbd solver (obj_). For DQ, obj_ and obj holds the value
    % before and after refinement, resp. For lasso, obj_ = obj.
    obj_trials = results.obj(idx,:);
    obj__trials = results.obj_(idx,:);
    its_trials = results.its(idx,:);

    if saveconvdata
      a0_trials = results.a0(idx,:);
      x0_trials = results.x0(idx,:);
      a_trials = results.a(idx,:);
    end

    % WHAT HAPPENS IN EACH TRIAL:
    t0 = tic;
    %for trial = 1:trials
    parfor (trial = 1:trials, n_workers)
      % A) Generate data: supp(x) must be >= 1
      a0 = adist(v1i, v2j);  %#ok<PFBNS>
      a0 = a0/norm(a0); 
      
      xgood = false;
      while ~xgood
        x0 = xdist(v1i, v2j); %#ok<PFBNS>
        xgood = sum(x0~=0) >= 1;
      end
      y = cconv(a0, x0, numel(x0));

      % B) Create solver and run continuation sequence
      solver = solverfun(y, a0, v1i, v2j); %#ok<PFBNS>
      [solver, stats] = solver.solve();

      % C) Record statistics
      obj_trials(trial) = maxdotshift(a0, solver.a);
      obj__trials(trial) = maxdotshift(a0, stats(end).a);
      its_trials(trial) = solver.it;

      if saveconvdata
        a0_trials{trial} = a0;
        x0_trials{trial} = sparse(x0);
        a_trials{trial} = solver.a;
      end
    end
    results.exp_idx = idx;
    results.times(idx) = toc(t0);

    % Write back to results
    results.obj(idx,:) = obj_trials;
    results.obj_(idx,:) = obj__trials;
    results.its(idx,:) = its_trials;

    if saveconvdata
      results.a0(idx,:) = a0_trials;
      results.x0(idx,:) = x0_trials;
      results.a(idx,:) = a_trials;
    end

    % Print & back up results
    fprintf('%s = %.2E, %s = %.2E, mean obj. = %.2f.', ...
      params.vars{2}, v2j, params.vars{1}, v1i, mean(results.obj(idx,:)));
    fprintf(' Time elapsed: %.1fs.\n', results.times(idx));
    if i == numel(theta);  disp(' ');  end

    if numel(backup);  assignin('base', backup, results);  end
  end

  % Final adjustments to results.
  results.obj = reshape(results.obj, [numel(theta) numel(p) trials]);
  results.obj_ = reshape(results.obj_, [numel(theta) numel(p) trials]);
  results.its = reshape(results.its, [numel(theta) numel(p) trials]);
  disp('Done.');
end

%#ok<*SAGROW>