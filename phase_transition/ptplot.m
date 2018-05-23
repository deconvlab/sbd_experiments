function ptplot(fname, pltcrit, lines, nticks)
load(fname, 'p0s', 'thetas', 'obj');

if nargin < 2 || isempty(pltcrit);  pltcrit = 0.9;      end
if nargin < 3 || isempty(nticks);   nticks = [6 6];     end

if isscalar(pltcrit)
    pltcrit = @(obj) mean(obj >= pltcrit, 3);
end

x = log10(thetas);  y = log10(p0s);

colormap gray;
imagesc(flipud(pltcrit(obj)'));
axis equal; 
xlim([0.5 size(obj,1)+0.5]);
ylim([0.5 size(obj,2)+0.5]);

xlabel('$\log_{10}(\theta)$','interpreter','latex','fontsize',16);
ylabel('$\log_{10}(p)$','interpreter','latex','fontsize',16);

set(gca,'TickLabelInterpreter','latex');
ticks = linspace(1,size(obj,1),nticks(1));
tlabels = linspace(x(1),x(end),nticks(1)); 
xticks(ticks);
xticklabels(tlabels);

ticks = linspace(1,size(obj,2),nticks(2));
tlabels = linspace(y(end),y(1),nticks(2));
yticks(ticks);
yticklabels(tlabels);

hold on;
lgd = cell(size(lines,1),1);
for l = 1:size(lines,1)
    xs = ([lines(l,1) lines(l,2)]-x(1))*(size(obj,1)+1)/(x(end)-x(1));
    ys = (y(end)-[lines(l,3) lines(l,4)])*(size(obj,2)+1)/(y(end)-y(1));
    plot(xs, ys, '--', 'linewidth', 3);
    
    lgd{l} = sprintf('$\\theta = p^{%.1f}$', ...
        (lines(l,2)-lines(l,1))/(lines(l,4)-lines(l,3)));
end
hold off;
legend(lgd, 'interpreter', 'latex', 'fontsize',14);
end

%#ok<*COLND>