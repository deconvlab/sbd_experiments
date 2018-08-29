classdef sbd_template < matlab.mixin.SetGet
properties
  a;  x;
  params = struct( ...
    'lambda', 0.1,  'alph', 0.9,  'data_init', true, ...
    'iter_lim', [1 10],  'iter_tol', 1e-3, ...
    'iter_lambdas', []
  );

  cost;
  it;
end

properties (Access = protected)
  y;  yhat;
  a0;  a_;
  s = obops;
end

methods
function o = sbd_template(y, ainit, params)
  if nargin >= 3 && ~isempty(params)
    o = set_params(o, params);
  end
  o = reset(set_y(o, y), ainit);
end

function o = set_y(o, y)
  o.y = y;
  o.yhat = fft(y);
end

function o = reset(o, ainit)
  m = numel(o.y);

  if nargin < 2 || isempty(ainit)
    ainit = numel(o.a0);
  elseif numel(ainit) > 1
    o.a0 = ainit;
  elseif o.params.data_init
    % MODIFY
    o.a0 = o.y(mod(randi(m) + (1:ainit), m) + 1);
    o.a0 = [zeros(ainit-1,1); o.a0(:); zeros(ainit-1,1)];
  else
    o.a0 = randn(ainit,1);
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
    costs = [Inf NaN(1, ilim(2))];
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

function o = center(o)
    p = numel(o.a);
    wsz = ceil((2+p)/3);

    [~, tau] = max(abs(xcorr(o.a, ones(wsz,1))));
    tau = tau - p - floor((p-wsz)/2);
    tmp = circshift([o.a; zeros(p,1)], -tau);
    o.a = tmp(1:p);
    o.x = circshift(o.x, tau);
end

% REPLACE WITH PARAMS
function [o, stats] = solve(o, ilim, tol, lambdas)
    o_lam = o.params.lambda;
    if nargin < 2 || isempty(ilim);     ilim = [];  end
    if nargin < 3 || isempty(tol);      tol = [];   end
    if nargin < 4 || isempty(lambdas)
        lambdas = o_lam;
    end

    stats = cell(1,numel(lambdas));
    for i = 1:numel(lambdas)
        o = center(o);
        o.params.lambda = lambdas(i);
        [o, stats{i}] = iterate(o, ilim, tol);
    end
    o.params.lambda = o_lam;
    stats = cell2mat(stats);
end

end
end