classdef sbd_dq < sbd_template

properties (Access = private)
  half_ynorm_sq;
  t
end

methods
function o = sbd_dq(y, ainit, params)
  o = o@sbd_template(y, ainit, params);
  o.half_ynorm_sq = norm(o.y)^2/2;
  o.t = 0.99/max(abs(o.yhat))^2;
end

function o = default_params(o)
  o.params = struct( ...
    'lambda', 0.1,  'alph', 0.9,  'data_init', true, ...
    'iter_ilim', [1 1e3],  'iter_tol', 1e-3, ...
    'solve_lambdas', [],  'solve_center', false, ...
    'refine_iters', 10*ones(10,2) ...
  );
end

function o = set_y(o, y)
  o = set_y@sbd_template(o, y);
  o.half_ynorm_sq = norm(o.y)^2/2;
end

function o = step(o)
  o.x = real(ifft(conj(o.yhat) .* fft(o.a, numel(o.y))));
  o.x = soft(o.x, o.params.lambda);
  xhat = fft(o.x);

  if o.params.alph > 0
      w = o.s.Exp(o.a, o.params.alph * o.s.Log(o.a_, o.a));
  else
      w = o.a;
  end
  g = -real(ifft(xhat .* o.yhat));
  g = g(1:numel(w));

  o.a_ = o.a;
  o.a = o.s.Exp(w, -o.t * o.s.e2rgrad(w, g));
  o.cost = o.half_ynorm_sq - norm(o.x)^2/2;
  o.it = o.it + 1;
end

function [o, stats] = solve(o)
  [o, stats] = solve@sbd_template(o);

  % Refinement
  function [agrad, ygrad] = refine_LS_grad(y,x,p) %#ok<DEFNU>
    ysupp = logical(cconv(x~=0, ones(numel(o.a),1), numel(x)));

    xhat = fft(x);
    ygrad = real(ifft(conj(xhat).*fft(ysupp.*y)));
    ygrad = ygrad(1:p);

    function out = agradfun(a)
      out = real(ifft(xhat.*fft(a, numel(x))));
      out = real(ifft(conj(xhat).*fft(ysupp.*out)));
      out = out(1:numel(a));
    end
    agrad = @agradfun;
  end

  n_iters = size(o.params.refine_iters,1);
  if n_iters
    lambda = 1;
    xsolver = set_y(lasso_fista(), o.y);
    xsolver.x = o.x;

    for i = 1:n_iters
      % Solve for x using lasso and get support
      xsolver = xsolver.set_params(struct('maxit', o.params.refine_iters(i,1)));
      xsolver = xsolver.evaluate(o.a, lambda);
      o.x = xsolver.x;

      % Update a using least squares / gradient descent on xsupp
      if false
        [agrad, ygrad] = refine_LS_grad(o.y, o.x, numel(o.a)); %#ok<UNRCH>
        o.a = pcg(agrad, ygrad);
      else
        a_ = o.a;
        ysupp = logical(cconv(o.x~=0, ones(numel(a_),1), numel(o.x)));
        xhat = fft(xsolver.x);
        L = max(abs(xhat))^2;
        for j = 1:o.params.refine_iters(i,2)
          agrad = o.a + 0.9 * (o.a - a_);
          agrad = real(ifft(xhat.*fft(agrad, numel(xhat))))-o.y;
          agrad = real(ifft(conj(xhat).*fft(ysupp.*agrad)));
          agrad = agrad(1:numel(a_));

          o.a = o.a - agrad/L;
        end
      end
      o.a = o.a/norm(o.a);

      % Reduce lambda
      lambda = lambda/2;
    end
  end
end

end
end