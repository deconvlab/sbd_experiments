classdef sbd_lasso < sbd_template
properties
    xsolver = lasso_fista(struct('maxit', 1e2, 'tol', 1e-4, 'xpos', false));
end

methods
function o = sbd_lasso(y, ainit, params, xparams)
  if nargin < 3 || isempty(params);  params = [];  end

  o = o@sbd_template(y, ainit, params);
  o.xsolver = set_y(o.xsolver, y);

  if nargin >= 4 && ~isempty(xparams)
      o.xsolver = set_params(o.xsolver, xparams);
  end
end

function o = set_y(o, y)
    o = set_y@sbd_template(o, y);
    o.xsolver = set_y(o.xsolver, y);
end

function o = reset(o, ainit)
    if nargin < 2;  ainit = [];  end
    o = reset@sbd_template(o, ainit);
    o.xsolver = reset(o.xsolver);
end

function o = step(o)
    o.xsolver.x = o.x;
    [o.xsolver, o.cost] = evaluate(o.xsolver, o.a, o.params.lambda);
    o.x = o.xsolver.x;
    xhat = fft(o.x);

    w = o.s.Exp(o.a, o.params.alph * o.s.Log(o.a_, o.a));
    g = real(ifft(conj(xhat) .* (xhat .* fft(w, numel(o.x)) - o.yhat)));
    g = g(1:numel(w));

    t = 0.99/max(abs(xhat))^2;
    o.a_ = o.a;
    o.a = o.s.Exp(w, -t*o.s.e2rgrad(w, g));
    o.it = o.it + 1;
end

end
end