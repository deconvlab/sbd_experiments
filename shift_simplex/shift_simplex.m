clc; clear; %#ok<*NOPTS>
run('../initpkg.m');

plotsurf = false;  plotcontour = true; %#ok<*NASGU>

%% Problem setup
p = 32;
a0 = randn(p,1);  a0 = a0/norm(a0);

m = 1024;                       % Observation size
theta = 2e-1;                   % Bernoulli (sparsity) coefficient
dist = @(m,n) randn(m,n);       % Distribution of activations
x0 = (rand(m,1) <= theta) .* dist(m,1);

y = cconv(a0, x0, m);

lambda = 0.1;

%% Calculate phi over the simplex
gsamps = linspace(-0.25,1.25, 1e2);        % coordinates on grid space

% Transformation from shift space to grid space -- x = Cu + d
C = [-.5 .5; -sqrt(3)/2 -sqrt(3)/2];
Cinv = inv(C);
d = [0.5; 0.5 + sqrt(3)/4];

% Shifts
s = 4;  s = -ceil(p/s):ceil(p/s);
s = s(randperm(numel(s), 3));
s = [0 s(s~=0)];  s = s(1:3);

u = [a0; zeros(m-p,1)];
v = circshift(u, s(3));  v = v(1:p);  v = v/norm(v) - a0;
u = circshift(u, s(2));  u = u(1:p);  u = u/norm(u) - a0;

phi_g = repmat({NaN(numel(gsamps),1)}, [1 numel(gsamps)]);
La = repmat({NaN(numel(gsamps),1)}, [1 numel(gsamps)]);
parfor i = 1:numel(gsamps)
    phi = set_y(lasso_fista(), y);
    for j = 1:numel(gsamps)
        uv = Cinv*([gsamps(i); gsamps(j)] - d); %#ok<MINV>
        
        a = a0 + uv(1)*u + uv(2)*v; 
        a = a/norm(a);
        [phi, phi_g{i}(j)] = evaluate(phi, a, lambda);
        La{i}(j) = max(abs(fft(phi.x)))^2;
    end
end
La = cell2mat(La);
phi_g = cell2mat(phi_g);
phidelta = 10 - min(phi_g(:));
phi_g = phi_g + phidelta;

% Plot phi surface
clf; hold off;
surfscript; 

%% Trajectories over hemisphere
maxit = 2e3;

% Fix initial point
uv = [1; 1]/3;
a = a0 + uv(1)*u + uv(2)*v; 
a = a/norm(a);

solvers = {
    sbd_lasso(y, a, struct('alph',0, 'stepsz', 0.99/max(La(:))))
    sbd_lasso(y, a, struct('alph',0.9, 'stepsz', 0.99/max(La(:))))
    };

phi = arrayfun(@(~) evaluate(set_y(lasso_fista(), y), a, lambda), ...
    1:numel(solvers), 'UniformOutput', false);

xypath = repmat({[C*uv + d NaN(2,maxit)]}, [numel(solvers) 1]);
uvpath = repmat({[uv NaN(2,maxit)]}, [numel(solvers) 1]);
costs = repmat({[phi{1}.costs(end) NaN(1, maxit)]}, [numel(solvers) 1]);
%debug = repmat({NaN(2, maxit)}, [numel(solvers) 1]);
parfor j = 1:numel(solvers)
    for i = 2:maxit+1
        solvers{j} = iterate(solvers{j});
        uv = -[u v -solvers{j}.a]\a0;  uv = uv(1:2);
        
        uvpath{j}(:,i) = uv;
        xypath{j}(:,i) = C*uv + d;
        
        a = [u v]*uv + a0;  a = a/norm(a);
        [phi{j}, costs{j}(i)] = evaluate(phi{j}, a, lambda);
        
        %debug{j}(:,i-1) = [norm(a) dot(a,solvers{j}.A)];
    end
end

%% Plotting
clf; hold off;

lgd = {'RGD', 'ARGD'};
colors = [1 .4 .3; 1 0 1; 0 0.5 0];
sym = {'s', 'd', '^'};
pit = ceil(sum(costs{1}>=costs{1}(end)+1e-1) * 0.99);
%pidxs = 1:10:pit+1;
pidxs = round(linspace(1, pit+1, 20));

h = [];
if plotsurf
for j = 1:numel(solvers)
    plot3(xypath{j}(1,1:pit+1), xypath{j}(2,1:pit+1), ...
        costs{j}(1:pit+1)+phidelta + 1 + 3*(j-1), ...
        'LineWidth', 1.5, 'Color', colors(j,:)); hold on;
    
    h = [h, plot3(...
        xypath{j}(1,pidxs), xypath{j}(2,pidxs), ...
        costs{j}(pidxs)+phidelta + 1 + 3*(j-1), ...
        sym{j}, 'MarkerSize', 10,...
        'LineWidth', 1, 'Color', colors(j,:) ...
        )];
end
end

if plotcontour
for j = 1:numel(solvers)
    plot3(xypath{j}(1,1:pit+1), xypath{j}(2,1:pit+1), ...
        zeros(pit+1,1) + 2*(j-1), ...
        'LineWidth', 1.5, 'Color', colors(j,:)); hold on;
    
    h = [h, plot3(xypath{j}(1,pidxs), xypath{j}(2,pidxs), ...
        zeros(size(pidxs)) + 2*(j-1), ...
        sym{j}, 'MarkerSize', 10,...
        'LineWidth', 1, 'Color', colors(j,:))]; %#ok<*AGROW,*UNRCH>
end
end
surfscript;

legend(h(1:numel(solvers)), lgd(1:numel(solvers)), 'Location', 'southwest');

%view(vw{:});
view(0,90);  
xlim([gsamps(1) gsamps(end)]);  
ylim([gsamps(1) gsamps(end)]);
