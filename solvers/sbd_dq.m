classdef sbd_dq < sbd_template
properties (Access = private)
    half_ynorm_sq;
    t
end
    
methods
function o = sbd_dq(y, ainit, params)
    o = o@sbd_template(y, ainit);
    o.half_ynorm_sq = norm(o.y)^2/2;
    o.t = 0.99/max(abs(o.yhat))^2;
    
    if nargin >= 3 && ~isempty(params)
        o = set_params(o, params);
    end
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

end
end