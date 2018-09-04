function plot2var(results, pltcrit, lines, fitline, nticks) %#ok<*INUSL>

if nargin < 2 || isempty(pltcrit);  pltcrit = [0.99 0.1]; end
if nargin < 3 || isempty(lines);    lines = [];           end
if nargin < 4 || isempty(fitline);  fitline = false;      end
if nargin < 5 || isempty(nticks);   nticks = [6 6];       end

if isnumeric(pltcrit)
  pltcrit = @(obj, obj_) mean( ...
    (obj >= pltcrit(1)) & (abs(obj_-obj) < pltcrit(2)), 3 ...
  );
end

args = {'vars', 'obj', 'obj_'};
for i = 1:numel(args)
  eval([args{i} '= results.' args{i} ';']);
end
x = vars{1,2}; %#ok<*IDISVAR,*USENS>
y = vars{2,2};

clf;
colormap gray;
M = pltcrit(obj, obj_);
imagesc(flipud(M'));
axis equal;
xlim([0.5 size(obj,1)+0.5]);
ylim([0.5 size(obj,2)+0.5]);


xlabel(sprintf('$%s$', vars{1,1}),'interpreter','latex','fontsize',16);
ylabel(sprintf('$%s$', vars{2,1}),'interpreter','latex','fontsize',16);

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

if ~isempty(lines)
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
    w0 = [1 lines(1,:)];
    w = logisticfit([X(:) Y(:)], M(:), 1e2*w0);
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








