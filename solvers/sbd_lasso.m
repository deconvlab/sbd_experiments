classdef sbd_lasso < sbd_template
properties
    xsolver = lasso_fista(struct('maxit', 1e2, 'tol', 1e-4, 'xpos', false));
end

methods
function o = sbd_lasso(params, xparams)
  if nargin < 1;  params = struct();   end
  if nargin < 2;  xparams = struct();  end

  o = o@sbd_template(params);
  o.xsolver.set_params(xparams);
end

function o = reset(o)
    o = reset@sbd_template(o);
    o.xsolver.reset();
end

function o = step(o)
  del_yvars = o.mk_yvars();
  
  w = o.s.Exp(o.a, o.params.alph * o.s.Log(o.a_, o.a));
  [g, o.cost, xhat] = o.calc_grad(w);

  t = 0.99/max(abs(xhat))^2;
  o.a_ = o.a;
  o.a = o.s.Exp(w, -t*o.s.e2rgrad(w, g));
  o.it = o.it + 1;
  o.yvars = ~del_yvars;
end

function [g, cost, xhat] = calc_grad(o, a)
  del_yvars = o.mk_yvars();
  
  o.xsolver.y = o.y;
  o.xsolver.x = o.x;
  cost = o.xsolver.evaluate(a, o.params.lambda);
  o.x = o.xsolver.x;
  
  xhat = fft(o.x);
  g = real(ifft(conj(xhat) .* (xhat .* fft(a, numel(o.x)) - o.yhat)));
  g = g(1:numel(a));
  
  o.yvars = ~del_yvars;
end

end
end