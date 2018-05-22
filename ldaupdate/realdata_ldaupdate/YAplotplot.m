clear; load('realdata_rec.mat');
figure(1); clf; 

imagesc(Y); axis off; axis equal;

p = get(gca, 'Position');
h = axes('Parent', gcf, 'Position', [p(1)+.45 p(2)+.5 p(3)-.55 p(4)-.55]);

imagesc(h, A); 
rectangle('Position', [0.5 0.5 size(A)], ...
    'EdgeColor', 'r', 'LineWidth', 3);
axis equal; axis off;