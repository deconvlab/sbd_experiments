function ptplot(fname, pltcrit, lines, nticks, fitline)
load(fname, 'p0s', 'thetas', 'obj', 'logt', 'logp');

if nargin < 2 || isempty(pltcrit);  pltcrit = 0.9;      end
if nargin < 4 || isempty(nticks);   nticks = [6 6];     end
if nargin < 5 || isempty(fitline);  fitline = false;    end

if isscalar(pltcrit)
    pltcrit = @(obj) mean(obj >= pltcrit, 3);
end

if exist('logt','var'); x=logt;  else x=log10(thetas); end
if exist('logp','var'); y=logp;  else y=log10(p0s);    end

clf;
colormap gray;
M = pltcrit(obj);
imagesc(flipud(M'));
axis equal;
xlim([0.5 size(obj,1)+0.5]);
ylim([0.5 size(obj,2)+0.5]);


xlabel('$\log_{10}(\theta)$','interpreter','latex','fontsize',16);
ylabel('$\log_{10}(p)$','interpreter','latex','fontsize',16);

ax = gca;
ax.TickLabelInterpreter ='latex';

ticks = linspace(1,size(obj,1),nticks(1));
tlabels = linspace(x(1),x(end),nticks(1));
ax.XTick = ticks;
ax.XTickLabels = tlabels;

ticks = linspace(1,size(obj,2),nticks(2));
tlabels = linspace(y(end),y(1),nticks(2));
ax.YTick = ticks;
ax.YTickLabels = tlabels;

if nargin >= 3 && ~isempty(lines)
  % If using a two-point system, convert to 2-affine system
  %   [x1 x1 y1 y2]  ==>  x + a*y + b = 0
  if size(lines,2) == 4
    xs = lines(:,1:2);  ys = lines(:,3:4);
    lines = -(xs(:,2)-xs(:,1))/(ys(:,2)-ys(:,1));
    lines = [lines -(xs(:,1)+lines.*ys(:,1))];
  end
  
  % If FITLINE==TRUE, initialize with first line to fit the halfway line
  if fitline
    [X, Y] = meshgrid(x,y);
    wsharp = 1e2;
    w0 = [1 lines(1,:)];
    M_ = 1./(1+exp( -wsharp*(w0(1)*X+w0(2)*Y+w0(3) )));
    wsign = sign(sum(sum((M_-0.5) .* (M-0.5))));
    w = logisticfit([X(:) Y(:)], M(:), wsign*wsharp*w0);
    lines(1,:) = w(2:3)/w(1);
  end
  
  % Calculate line in p-theta space
  xs = (min(x) + max(x))/2;
  xs = ([min(x) max(x)] - xs)*2 + xs;
  ys = -(xs + repmat(lines(:,2), [1 2]))./repmat(lines(:,1), [1 2]);
  
  % Convert line to grid space
  %   [x y]' = [x0 y0]' + [cu dv]'
  %   [c d]' = [dx dy]' ./ [du dv]'
  us = (xs - min(x))*(numel(x)-1)/range(x)+1;
  vs = -(ys - max(y))*(numel(y)-1)/range(y)+1;
  
  hold on;
  plot(us, vs', '--', 'linewidth', 3);
  lgd = arrayfun(@(c) sprintf('$\\theta \\simeq p^{%.3f}$',c), ...
    -lines(:,1), 'UniformOutput', false);
  legend(lgd, 'interpreter', 'latex', 'fontsize',14);
  hold off;
end
end

%#ok<*COLND>








