classdef sbd_dq < sbd_template

properties (Access = protected)
  half_ynorm_sq;
  tmin;
end

methods
function o = sbd_dq(params)
  if nargin < 1;  params = [];  end
  o = o@sbd_template(params);
end

function o = default_params(o)
  o.params = struct( ...
    'lambda', 0.1,  'alph', 0.9, ...
    'iter_ilim', [1 1e3],  'iter_tol', 1e-3, ...
    'solve_lambdas', [],  'solve_center', false, ...
    'backtrack', [0.1 0.1], ...
    'refine_iters', 10*ones(10,2) ...
  );
end

function del_tmpvars = mk_tmpvars(o)
% Only create and delete tmpvars at the outermost loop.
  del_tmpvars = ~o.tmpvars;    
  if ~o.tmpvars   % only happens if tmpvars are invalid
    mk_tmpvars@sbd_template(o);
    o.half_ynorm_sq = norm(o.y)^2/2;
    o.tmin = 0.99/max(abs(o.yhat))^2;
  end
end

function [g, x] = calc_grad(o, a)
  del_tmpvars = o.mk_tmpvars();  
  
  x = real(ifft(o.yhat .* conj(fft(a, numel(o.y)))));
  x = soft(x, o.params.lambda);
  xhat = fft(x);

  g = -real(ifft(conj(xhat) .* o.yhat));
  g = g(1:numel(a));
  
  o.tmpvars = ~del_tmpvars;
end

function o = step(o)
  del_tmpvars = o.mk_tmpvars();  
  
  if o.params.alph > 0
      w = o.s.Exp(o.a, o.params.alph * o.s.Log(o.a_, o.a));
  else
      w = o.a;
  end
  [g, x] = o.calc_grad(w);
  g = o.s.e2rgrad(w, g);
  o.gnorm2 = norm(g)^2;
  cost_ = o.half_ynorm_sq - norm(x)^2/2;

  if ~isempty(o.params.backtrack)
    o.t = 1;  btdec = o.params.backtrack(1);  btslack = o.params.backtrack(2);
  else
    o.t = o.tmin;  btslack = 1;
  end

  repeat = true;
  while repeat
    o.t = max(o.t, o.tmin); %#ok<*PROP>
    a = o.s.Exp(w, -o.t * g);
    [~, x] = o.calc_grad(a);
    cost = o.half_ynorm_sq - norm(x)^2/2;

    repeat = (cost-cost_ >= -btslack*o.t*o.gnorm2) && (o.t > o.tmin);
    if repeat;  o.t = btdec * o.t;  end  %else;  disp(o.t/o.tmin);  end
  end

  o.a_ = o.a;  o.a = a;  o.x = x;
  o.cost = cost;
  o.it = o.it + 1;
  o.tmpvars = ~del_tmpvars;
end

function [o, stats] = solve(o)
  del_tmpvars = o.mk_tmpvars();  
  
  % Iterate for all lambdas
  [o, stats] = solve@sbd_template(o);

  % Refinement
  n_refine = size(o.params.refine_iters,1);
  if n_refine
    lambda = 1;
    xsolver = set_y(lasso_fista(), o.y);
    xsolver.x = o.x;

    for i = 1:n_refine
      % Solve for x using lasso and get support
      xsolver = xsolver.set_params(struct('maxit', o.params.refine_iters(i,1)));
      xsolver = xsolver.evaluate(o.a, lambda);
      o.x = xsolver.x;

      % Update a using gradient descent on xsupp
      a_ = o.a;
      ysupp = logical(cconv(o.x~=0, ones(numel(a_),1), numel(o.x)));
      xhat = fft(xsolver.x);
      L = max(abs(xhat))^2;
      for j = 1:o.params.refine_iters(i,2)
        agrad = o.a + 0.9 * (o.a - a_);
        agrad = real(ifft(xhat.*fft(agrad, numel(xhat))))-o.y;
        agrad = real(ifft(conj(xhat).*fft(ysupp.*agrad)));
        agrad = agrad(1:numel(a_));

        a_ = o.a;
        o.a = o.a - agrad/L;
      end
      o.a = o.a/norm(o.a);

      % Reduce lambda
      lambda = lambda/2;
    end
  end
  o.tmpvars = ~del_tmpvars;
end

end
end