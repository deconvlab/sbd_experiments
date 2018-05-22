tmp = {0.95 'flat'};
i = log10(thetas);  j = log10(p0s);         % *params

clf;
if min(size(obj,1), size(obj,2)) == 1
    % 1D transition plots
    if size(obj,1) == 1;  tmp{3} = j;  else;  tmp{3} = i;  end
    yyaxis left;
    plot(tmp{3}, median(squeeze(obj),2));  hold on;
    plot(tmp{3}, mean(squeeze(obj),2));
    plot(tmp{3}, min(squeeze(obj), [], 2));
    hold off;
    ylabel('\rho = max_i <S_i[a0], a>');

    yyaxis right;    plot(tmp{3}, mean(squeeze(obj) >= tmp{1},2));
    ylabel(['P[\rho > ' sprintf('%.2f]', tmp{1})]);
    xlabel('p\theta');

else
    % 2D transition plots
    colormap gray;
    subplot(131); surf(i, j, mean(obj,3)');
    view(2); shading(tmp{2});  title('Mean of \rho');
    xlim([i(1) i(end)]); ylim([j(1) j(end)]);
    xlabel('log(\theta)'); ylabel('log(p)');

    subplot(132); surf(i, j, median(obj,3)');
    view(2); shading(tmp{2});  title('Median of \rho');
    xlim([i(1) i(end)]); ylim([j(1) j(end)]);
    xlabel('log(\theta)'); ylabel('log(p)');

    subplot(133); surf(i, j, mean(obj>=tmp{1},3)');
    view(2); shading(tmp{2});  title('Success prob.');
    xlim([i(1) i(end)]); ylim([j(1) j(end)]);
    xlabel('log(\theta)'); ylabel('log(p)');
end