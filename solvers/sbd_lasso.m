classdef sbd_lasso < sbd_template
properties
    xsolver = lasso_fista(struct('maxit', 1e2, 'tol', 1e-4));
end

methods
function o = sbd_lasso(y, ainit, params, xparams)
    o = o@sbd_template(y, ainit);
    reset(o);
    
    if nargin >= 3 && ~isempty(params)
        o = set_params(o, params);
    end
    
    if nargin >= 4 && ~isempty(xparams)
        o.xsolver = set_params(o.xsolver, xparams);
    end
end

function o = reset(o, y, ainit)
    if nargin < 2;  y = [];  end
    if nargin < 3
        reset@sbd_template(o, y);
    else
        reset@sbd_template(o, y, ainit);
    end
    o.xsolver = reset(o.xsolver);
end

function o = step(o)
    [o.xsolver, o.cost] = evaluate(o.xsolver, o.y, o.a, o.params.weights);
    o.x = o.xsolver.x;
    xhat = fft(o.x);
    
    w = o.s.Exp(o.a, o.params.alph * o.s.Log(o.a_, o.a));
    g = real(ifft(conj(xhat) .* (xhat .* fft(w, numel(o.y)) - fft(o.y))));
    g = g(1:numel(w));
    
    t = 0.99/max(abs(xhat))^2;
    o.a_ = o.a;
    o.a = o.s.Exp(w, -t*o.s.e2rgrad(w, g));
    o.it = o.it + 1;
end

end
end