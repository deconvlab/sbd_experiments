%#ok<*SAGROW>
warning('OFF', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
for idx = idx0:prod(tmp)-1
  fprintf('[%5d|%5d]:  ', idx, prod(tmp)-1);
  idx_1 = idx+1;
  [i, j] = ind2sub(tmp, idx_1);

  theta = thetas(i); p0 = p0s(j);         % *params
  a0 = randn(p0,1);  a0 = a0/norm(a0);

  m = 1e2 * p0;
  lambda = [1 0.9] * 0.8/sqrt(p0*theta);

  start = tic;

  % WHAT HAPPENS IN EACH TRIAL:
  for trial = 1:trials
    % A) Generate x & y: supp(x) must be >= 1
    xgood = false;
    while ~xgood
      x0 = (rand(m,1) <= theta) .* dist(m,1);
      xgood = sum(x0~=0) >= 1;
    end
    y = cconv(a0, x0, m);

    % B) Create solver and run continuation sequence
    solver = solverfun(y, p0, struct('lambda', lambda));
    solver.solve([10 maxit], tol, lambda);

    % C) Record statistics
    obj(idx_1, trial) = maxdotshift(a0, solver.a);
    its(idx_1, trial) = solver.it;
  end

  fprintf('p0 = %d, theta = %.2E, mean obj. = %.2f.', ...
      p0, theta, mean(obj(idx_1,:)));       % *params
  times(idx_1) = toc(start);
  fprintf(' Time elapsed: %.1fs.\n', times(idx_1));
  if i == numel(thetas);  disp(' ');  end
  % save('tmp.mat');
end

%%
warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
obj = reshape(obj, [tmp trials]);
its = reshape(its, [tmp trials]);
disp('Done.');