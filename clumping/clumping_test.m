clear; clc;
run('../initpkg.m');

%% Data
n = 5e2;
p = 10;

theta = 0.3;
theta0 = 0.1;
r1 = 1;

x0 = sqclumpx(n, [0.1 (theta-(1-r1)*theta0)/r1 r1 5]);
a0 = randn([p 1]);
a0 = a0/norm(a0);

y = cconv(x0, a0, n);

%% Deconvolve test
solver = sbd_lasso(y, p, struct('iter_ilim', [1 1000]), struct('xpos', true));
solver.iterate();

subplot(1,3,2:3); stem([x0 solver.x], '.')
subplot(131); stem([[a0 ; zeros(numel(solver.a)-p,1)] solver.a], '.');