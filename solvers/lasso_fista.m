classdef lasso_fista < handle
properties
    x;  y;
    costs;
    params = struct('maxit', 1e4, 'tol', 1e-4, 'xpos', false);

    it;
    eps;
end

methods
function o = lasso_fista(params)
  if nargin < 1;  params = struct();  end  
  o = set_params(o, params);
end

function o = reset(o)
    o.x = [];
    o.costs = [];
    o.it = [];
    o.eps = [];
end

function o = set_params(o, params)
    tmp = intersect(fieldnames(o.params), fieldnames(params));
    for i = 1:numel(tmp)
        o.params.(tmp{i}) = params.(tmp{i});
    end
end

function cost = evaluate(o, a, weights)
    m = numel(o.y);  a = a(:);

    ahat = fft(a,m);
    a2hat = abs(ahat).^2;
    ayhat = conj(ahat).*fft(o.y);
    s = 0.99/max(a2hat);

    if isempty(o.x)
        o.x = randn(m,1);
    end

    w = fft(o.x);  xhat_ = w;  t = 1;
    o.it = 0;  repeat = true;
    o.costs = NaN(o.params.maxit,1); cost = Inf;
    while repeat
        o.x = soft(real(ifft(w - s*(a2hat.*w-ayhat))), s*weights);
        if o.params.xpos;  o.x = max(o.x, 0);  end
        
        t_ = (1+sqrt(1 + 4*t^2))/2;
        xhat = fft(o.x);
        w = xhat + (t-1)/t_ * (xhat - xhat_);

        t = t_;  xhat_ = xhat;
        o.it = o.it + 1;

        o.costs(o.it) = norm(real(ifft(ahat.*xhat))-o.y)^2/2 ...
            + norm(weights .* o.x(:),1);

        o.eps = abs(cost - o.costs(o.it));
        cost = o.costs(o.it);
        repeat = (o.it < o.params.maxit) && (o.eps >= o.params.tol);
    end

    o.costs = o.costs(1:o.it);
end
end
end