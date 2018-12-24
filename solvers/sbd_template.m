classdef sbd_template < matlab.mixin.SetGet
properties
  a;  x;  y;
  params;

  cost;
  it;
end

properties (Access = protected)
  p0;
  yvars = false;  
  yhat;
  ainit;  a_;
  s = obops;
end

methods
function o = sbd_template(params)
  if nargin < 1;  params = struct();  end
  o = o.default_params();
  o.set_params(params);
  o.reset();
end

function o = default_params(o)
  o.params = struct( ...
    'lambda', 0.1,  'alph', 0.9, ...
    'iter_ilim', [1 1e3],  'iter_tol', 1e-3, ...
    'solve_lambdas', [],  'solve_center', false ...
  );
end

function o = set_params(o, params)
  tmp = intersect(fieldnames(o.params), fieldnames(params));
  for i = 1:numel(tmp)
    o.params.(tmp{i}) = params.(tmp{i});
  end
end

function o = reset(o)
  o.a = o.ainit;  o.a_ = o.ainit;
  o.x = [];
  o.it = 0;
  o.cost = [];
end

function o = set_ainit(o, ainit)
  if nargin < 2 || isempty(ainit)
    ainit = o.p0;
  end

  if numel(ainit) > 1
    o.ainit = ainit;
    o.p0 = numel(ainit);
    o.ainit = o.ainit(:)/norm(o.ainit(:));
  else
    o.ainit = data_init(o, ainit);
  end

  o.a = o.ainit;  o.a_ = o.ainit;
end

function ainit = data_init(o, p0, calc_grad)
  if nargin < 3;  calc_grad = false;  end
  m = numel(o.y);

  ainit = o.y(mod(randi(m) + (1:p0), m) + 1);
  ainit = [zeros(p0-1,1); ainit(:); zeros(p0-1,1)];
  if calc_grad
    ainit = -o.calc_grad(ainit/norm(ainit(:)));
  end
  ainit = ainit(:)/norm(ainit(:));
end

function del_yvars = mk_yvars(o)
% Only create and delete yvars at the outermost loop.
  del_yvars = ~o.yvars;    
  if ~o.yvars   % only happens if yvars are invalid
    o.yvars = true;
    o.yhat = fft(o.y);
  end
end

function o = step(o)
end

function g = calc_grad(o, a) %#ok<INUSD,STOUT>
end

function [o, stats] = iterate(o)
  del_yvars = o.mk_yvars();
  assert(~isempty(o.y), 'y has not been initialized.');
  assert(~isempty(o.a), 'a has not been initialized.');
  
  ilim = o.params.iter_ilim;
  tol = o.params.iter_tol;

  repeat = true;
  i = 0;
  costs = [Inf NaN(1, ilim(2))];
  while repeat
    o = step(o);
    i = i + 1;
    costs(i+1) = o.cost;

    eps = abs(costs(i+1) - costs(i));
    repeat = (i < ilim(1)) || ((i < ilim(2)) && (eps > tol));
  end
  stats.a = o.a;
  stats.it = i;
  stats.eps = eps;
  stats.costs = costs(2:i+1);
  o.yvars = ~del_yvars;
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

function [o, stats] = solve(o)
  del_yvars = o.mk_yvars();
  o_lam = o.params.lambda;

  if ~isempty(o.params.solve_lambdas)
    lambdas = o.params.solve_lambdas;
  else
    lambdas = o_lam;
  end

  stats = cell(1,numel(lambdas));
  for i = 1:numel(lambdas)
    if o.params.solve_center
      o = center(o);
    end

    o.params.lambda = lambdas(i);
    [o, stats{i}] = iterate(o);
  end
  o.params.lambda = o_lam;
  stats = cell2mat(stats);
  o.yvars = ~del_yvars;
end

end
end