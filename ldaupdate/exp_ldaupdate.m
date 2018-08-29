%% Experiment Parameters
run('initpkg.m');
testcase = 'dq';
nexp = 20;
tol = 1e-20;
dist = @(m,n) randn(m,n); 

if strcmp(testcase, 'lasso')
    p = 5e2;
    m = p*1e2;
    theta = p^(-1/2);
    maxit = 200;
    eta = 1.01;
end

if strcmp(testcase, 'dq')
    p = 5e3;
    m = p*1e2;
    theta = p^(-2/3);
    maxit = 800;
    eta = 1.005;
end


%% Compare convergence rate with different 
for I = 1:nexp
    a0 = randn(p,1);
    a0 = a0/norm(a0);
    x0 = (rand(m,1) <= theta) .* dist(m,1);
    y = cconv(a0, x0, m);
    lambda = 0.8/sqrt(p*theta);
    
    %-change lambda
    if strcmp(testcase,'lasso')
        solver = sbd_lasso(y, 3*p-2, struct('alph',0,'lambda',lambda));
    elseif strcmp(testcase,'dq')
        solver = sbd_dq(y, 3*p-2, struct('alph',0,'lambda',lambda));
    end
    ainit = solver.a;
    for J = 1:maxit
        solver.params.lambda = min(solver.params.lambda*eta,0.5);
        solver = step(solver);
        if strcmp(testcase,'lasso')
            stats_lasso_var.dists(I,J) = maxdotshift(a0,solver.a,0);
        elseif strcmp(testcase,'dq')
            stats_dq_var.dists(I,J) = maxdotshift(a0,solver.a,0);
        end
    end
    
    %-fix lambda
    solver = reset(solver, ainit);
    for J = 1:maxit
        solver.params.lambda = lambda;
        solver = step(solver);
        if strcmp(testcase,'lasso')
            stats_lasso_fix.dists(I,J) = maxdotshift(a0,solver.a,0);
        elseif strcmp(testcase,'dq')
            stats_dq_fix.dists(I,J) = maxdotshift(a0,solver.a,0);
        end
    end
end

save exp_ldaupdate stats_lasso_var stats_lasso_fix ...
                   stats_dq_var    stats_dq_fix     nexp

%% Plot Result 
figure(1); clf; hold on; box on;
plot(sum(2-2*stats_lasso_fix.dists,1)/20,'x','linewidth',3);
plot(sum(2-2*stats_lasso_var.dists,1)/20,'x','linewidth',3);
plot(sum(2-2*stats_dq_fix.dists,1)/nexp,'linewidth',3);
plot(sum(2-2*stats_dq_var.dists,1)/nexp,'linewidth',3);


set(gca,'TickLabelInterpreter','latex', ...
        'linewidth',2, ...
        'fontsize',16,...
        'xtick',[0,200,500,800],...
        'ytick',[0,0.1,0.4,0.7,1.0,1.3],...
        'yticklabel',{'0.0',0.1,0.4,0.7,'1.0',1.3},...
        'ylim',[0,1.5]);
legend({'$(P_{\mathrm{lasso}})\lambda_{(i)} = \lambda_{(0)}$','$(P_{\mathrm{lasso}}),\lambda_{(i)} = 1.01\lambda_{(i-1)}$'...
         '$(P_{\mathrm{nc}}),\lambda_{(i)} = \lambda_{(0)}$','$(P_{\mathrm{nc}}),\lambda_{(i)} = 1.005\lambda_{(i-1)}$' },...
        'interpreter','latex','location','northeast','fontsize',14 );
xlabel('iteration number','interpreter','latex','fontsize',16);
ylabel('$d^2(a,\mathcal{D})$','interpreter','latex','fontsize',16);
plot(1:800,0.1*ones(800,1),'--k');
hold off;

%% Done
beep
