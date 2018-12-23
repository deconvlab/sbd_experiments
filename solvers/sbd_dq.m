classdef sbd_dq < sbd_template

properties (Access = private)
  half_ynorm_sq;
  t
end

methods
function o = sbd_dq(y, params)
  o = o@sbd_template(y, params);
  o = set_y(o,y);
end

function o = default_params(o)
  o.params = struct( ...
    'lambda', 0.1,  'alph', 0.9, ...
    'iter_ilim', [1 1e3],  'iter_tol', 1e-3, ...
    'solve_lambdas', [],  'solve_center', false, ...
    'backtrack', NaN, ...
    'refine_iters', 10*ones(10,2) ...
  );
end

function o = set_y(o, y)
  o.y = y;
  o.yhat = fft(y);
  o.half_ynorm_sq = norm(o.y)^2/2;
  o.t = 0.99/max(abs(o.yhat))^2;
end

function [g, x] = calc_grad(o, a)
  x = real(ifft(o.yhat .* conj(fft(a, numel(o.y)))));
  x = soft(x, o.params.lambda);
  xhat = fft(x);

  g = -real(ifft(conj(xhat) .* o.yhat));
  g = g(1:numel(a));
end

function o = step(o)
  if o.params.alph > 0
      w = o.s.Exp(o.a, o.params.alph * o.s.Log(o.a_, o.a));
  else
      w = o.a;
  end
  [g, x] = o.calc_grad(w);
  cost_ = o.half_ynorm_sq - norm(x)^2/2;

  o.a_ = o.a;
  bt = o.params.backtrack;
  if (0 < bt) && (bt < 1);  t = 1;  else;  t = o.t;  end

  repeat = true;
  while repeat
    t = max(t, o.t); %#ok<*PROP>
    o.a = o.s.Exp(w, -t * o.s.e2rgrad(w, g));
    [~, o.x] = o.calc_grad(w);
    cost = o.half_ynorm_sq - norm(o.x)^2/2;

    repeat = (cost-cost_ >= -t*norm(g)) && (t > o.t);
    if repeat;  t = bt * t;  end
  end

  o.cost = cost;
  o.it = o.it + 1;
end

function [o, stats] = solve(o)
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
end

end
end