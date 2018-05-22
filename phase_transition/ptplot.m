function ptplot(fname, pltcrit, nticks)
load(fname, 'p0s', 'thetas', 'obj');

if nargin <= 2 || isempty(pltcrit); pltcrit = 0.9;  end
if nargin <= 3 || isempty(nticks);  nticks = 6;     end

if isscalar(pltcrit)
    pltcrit = @(obj) mean(obj >= pltcrit, 3);
end

colormap gray;
imagesc(flipud(pltcrit(obj)'));
axis equal; 
xlim([0.5 size(obj,1)+0.5]);
ylim([0.5 size(obj,2)+0.5]);

xlabel('$\log_{10}(\theta)$','interpreter','latex','fontsize',16);
ylabel('$\log_{10}(p)$','interpreter','latex','fontsize',16);

set(gca,'TickLabelInterpreter','latex');
ticks = linspace(1,size(obj,1),nticks);
tlabels = linspace(log10(thetas(1)),log10(thetas(end)),nticks); 
xticks(ticks);
xticklabels(tlabels);

ticks = linspace(1,size(obj,2),nticks);
tlabels = linspace(log10(p0s(end)),log10(p0s(1)),nticks);
yticks(ticks);
yticklabels(tlabels);

hold on;
r = log10(p0s(end)/p0s(1))/log10(thetas(end)/thetas(1));
plot([0 size(obj,1)+1], [0 3/4 * r*size(obj,2)+1], 'r--', 'linewidth', 3);
plot([0 size(obj,1)+1], [0 1/2 * r*size(obj,2)+1], '-.', 'linewidth', 3);
hold off;
legend({'$\theta = p^{-3/4}$', '$\theta = p^{-1/2}$'}, ...
    'interpreter', 'latex', 'fontsize',14);
end

%#ok<*COLND>