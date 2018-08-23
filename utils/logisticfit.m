function [w, costs, yhat, info] = logisticfit(X, y, w0, params)
%LOGISTICFIT    Fits a logistic function to a single response variable.
%   Fits a logistic function to minimize the squared error of the
%   reconstruction for a single response variable. Designed for finding
%   the halfway point for affine transition regions in phase plots.
%
%   Because the gradients of the loss function don't admit analytic
%   roots or simple smoothness calculations, multiple local minimizers are
%   found using accelerated gradient descent with backtracking. The
%   coefficients producing the lowest cost is returned.
%
%   [W, COSTS, YHAT, INFO] = LOGISTICFIT(X, Y, W0, PARAMS)
%
%   Inputs:
%     X: (N x P) float. P predictors for each of the N data points.
%
%     Y: (N) float.  Response variable for each data point.
%
%   Optional inputs:
%     W0: (P+1) float.  
%       Initializations for the logistic coefficients. The final entry 
%       gives the bias s.t. the halfway point of the logistic model is 
%       given by X*W0(1:end-1) + W0(end). By default W0 = [1 .... 1 0]'.
%
%       The cost function contains many local minima so a good
%       initialization is important. The halfway point is constant for any
%       nonzero scaling of W0, but larger values for W0 correspond to a
%       sharper transition and often leads to better estimates.
%   
%     PARAMS: struct.  a number of hyperparameters for the function:
%       'alph': float [0,1).  Momentum parameter. 0.9 by default.  
%       
%       'btparams': (2) float.  The backtracking decay rate and minimum
%       stepsize, respectively.  [1e-1 1e-50] by default.
%       
%       'getb': bool.  Whether to estimate the bias. The bias is not
%       updated if 'getb' == FALSE. Typically this leads to better
%       estimates of the other coefficients, e.g. the slope.
%       
%       'niter': int.  Number of iterations to run for.
%
%       'trials': [int float].  Number of trials to take and additive noise 
%       std to W0 for each trial, respectively. [100 5] by default.       
%
%       't0': float.  Initial step size.
%
%
%   Returns:
%     W: (P+1) float.  The logistic model coefficients corresponding to the
%       minimum cost trial ; see W0.
%
%     COSTS: (NITER+1 x NTRIALS) float.  The costs of each trial over the 
%       course of optimization. NITER(1,:) gives the initial cost with 
%       initialized coefficients.
%
%     YHAT: (N) float.  Reconstructed response corresponding to the minimum
%       cost trial.
%
%     INFO: struct.  Contains W and YHAT across *all* trials and
%       backtracking information.
%
%  Written by Yenson Lau: 
%   https://github.com/yenson-lau/matlab_logisticfit.git
%

% Process prameters
niter = 500;
getb = false;
trials = [100 5];
alph = 0.99;
btparams = [1e-1 1e-50];
t0 = 1;
testgrad = false;
if nargin >= 4 && isstruct(params)
  paramnames = {'niter', 'getb', 'trials', 'alph', 'btparams', 't0', 'testgrad'};
  p = fieldnames(params);
  for i = 1:numel(p)
    assert(ismember(p{i}, paramnames), 'Invalid parameter specified.');
    eval([p{i} '= params.(p{i});']);
  end
end

% Initialize coefficients and bias over several trials with noise
if nargin < 3 || isempty(w0)
  w = ones(size(X,2),1);
  b = 0;
else
  w = w0(1:end-1);  w = w(:);
  b = w0(end);
end
w = repmat(w, [1 trials(1)]) + trials(2)*randn(size(w,1), trials(1));
b = repmat(b, [1 trials(1)]) + trials(2)*randn(1, trials(1));

% Detect the correct sign for the coefficients
costs = zeros(niter+1,trials(1));
[costs(1,:), r, yhat, EXP] = cost(X,y,w,b);
j = sign(sum(sum( (yhat-0.5).*repmat(y-0.5,[1, trials(1)]) )));
if j < 0
    w = -w;  b = -b;
    [costs(1,:), r, yhat, EXP] = cost(X,y,w,b);
end

% Initialize loop and iterate
t = t0 * ones(1,trials(1));
smlsteps = zeros(niter,trials(1));
w1 = w;  b1 = b;
gb = zeros(1,trials(1)); %#ok<PREALL>
for i = 1:niter
  % Calculate gradients
  J = r.*EXP./(1+EXP).^2;
  gw = X'*J;
  if getb;  gb = sum(J,1);  end;

  if testgrad
    j = randi(trials(1)); %#ok<UNRCH>
    testgrads(X,y,w1(:,j),b1(j),gw(:,j),gb(j),3,1e-1);
  end

  % Backtracking
  w_ = w;  b_ = b;
  bt = true(1, trials(1));
  cost_ = zeros(1, trials(1));
  while sum(bt)
    % Gradient step
    w(:,bt) = w1(:,bt) - repmat(t(bt), [size(w,1) 1]).*gw(:,bt);
    if getb;  b(bt) = b1(bt) - t(bt).*gb(bt);  end

    % Check costs and decrease stepsizes with large cost
    cost_(bt) = cost(X,y,w(:,bt),b(bt));
    lrgcost = cost_ > costs(i,:);
    t(lrgcost) = t(lrgcost) * btparams(1);

    % Stop backtracking when stepsize is too small or cost has decreased
    smlstep = t < btparams(2);
    bt(~lrgcost|smlstep) = false;
  end
  % Update any variables that decreased cost, keep rest the same
  costs(i+1,:) = cost_;
  w(:,smlstep) = w_(:,smlstep);
  b(:,smlstep) = b_(:,smlstep);
  costs(i+1,smlstep) = costs(i,smlstep);
  smlsteps(i,:) = smlstep;

  % Add momentum for next iteration; compute cost
  w1 = w + alph * (w - w_);
  if getb;  b1 = b + alph * (b - b_);  end
  [~, r(:,bt), yhat(:,bt), EXP(:,bt)] = cost(X,y,w1(:,bt),b1(bt));
end

% Save outputs
info.W = [w; b];
info.Yhat = yhat;
info.smlsteps = smlsteps;

[~,i] = min(costs(end,:));
w = info.W(:,i)';
yhat = info.Yhat(:,i);

end

function [c, r, yhat, EXP] = cost(X, y, w, b)
%COST    Calculate cost / intermediate variable
EXP = exp(-(X*w + repmat(b, [size(X,1) 1])));   % main expense: exp()
yhat = 1./(1+EXP);
r = yhat - repmat(y, [1 size(w,2)]); 
c = sum(r.^2,1)/2;
end

function testgrads(X, y, w, b, gw, gb, nw, delta)  %#ok<*DEFNU>
%TESTGRADS    Function for testing gradients
cost0 = cost(w,b);

N = 1000;
delta = linspace(-delta,delta,N)';

% Unit vectors
dw = randn(numel(w), nw);
for i = 1:nw
    dw(:,i) = dw(:,i)/norm(dw(:,i));
end

% Linearized cost estimates
dfw = gw'*dw;
linw = cost0 + delta*dfw;
linb = cost0 + delta*gb;

% Costs
costsw = zeros(N,nw);
costsb = zeros(N,1);
for i = 1:N
  for j = 1:nw
    costsw(i,j,1) = cost(X, y, w+delta(i)*dw(:,j), b);
  end
  costsb(i) = cost(X, y, w, b+delta(i));
end

subplot(121);
plot(delta,costsw); hold on;
set(gca,'ColorOrderIndex',1)
plot(delta, linw, '--'); hold off;

subplot(122);
plot(delta, costsb); hold on;
set(gca,'ColorOrderIndex',1)
plot(delta, linb, '--'); hold off;
pause;
end
