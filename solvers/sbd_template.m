classdef sbd_template < matlab.mixin.SetGet
properties
    y;  a;  x;
    params = struct('weights', 0.1, 'alph', 0.9);

    cost;
    it;
end

properties (Access = protected)
    a0;  a_;
    s = obops;
end

methods
function o = sbd_template(y, ainit)
    reset(o, y, ainit);
end

function o = reset(o, y, ainit)
    if nargin >= 2 && ~isempty(y)
        o.y = y;
    end
    if nargin < 3 || isempty(ainit)
        ainit = numel(o.a0);
    end
    
    if numel(ainit) > 1
        o.a0 = ainit;
    else
        tmp = mod(randi(numel(o.y)) + (1:ainit), numel(o.y)) + 1;
        o.a0 = o.y(tmp);
    end
    
    o.a0 = o.a0(:)/norm(o.a0(:));
    o.a = o.a0;  o.a_ = o.a0;
    o.x = [];
    
    o.it = 0;
    o.cost = [];
end

function o = set_params(o, params)
    tmp = intersect(fieldnames(o.params), fieldnames(params));
    for i = 1:numel(tmp)
        o.params.(tmp{i}) = params.(tmp{i});
    end
end

function o = step(o)
end

function [o, stats] = iterate(o, ilim, tol)
    if nargin < 2 || isempty(ilim);     ilim = [1 1];   end
    if nargin < 3 || isempty(tol);      tol = 0;        end
    
    repeat = true;
    it = 0;
    costs = [Inf; NaN(ilim(2),1)];
    while repeat
        o = step(o);
        it = it + 1;
        costs(it+1) = o.cost;
        
        eps = abs(costs(it+1) - costs(it));
        repeat = (it < ilim(1)) || ((it < ilim(2)) && (eps > tol));
    end
    stats.it = it;
    stats.eps = eps;
    stats.costs = costs(2:it+1);
end
end
end