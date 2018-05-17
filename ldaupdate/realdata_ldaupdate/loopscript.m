for i = i0:maxit
    stime = tic;
    parfor t = 1:trials
        s1{t} = iterate(s1{t});
        s2{t} = iterate(s2{t});
        
        s2{t}.f{1}.weights = min(s2{t}.f{1}.weights * eta, ldamax);
        s2{t}.f{1} = compile_params(s2{t}.f{1});
        
        A1s{t,i} = s1{t}.A{1};
        A2s{t,i} = s2{t}.A{1};
    end
    
    if ismember(i, updates)
        plotscript;
        
        fprintf(['Iter. %d.  '...
            'Costs: [%.2E %.2E]. Elapsed time %.2fs.\n'], ...
            i, s2{1}.cost, s2{1}.cost, toc(stime));
    end
    save('./realdata_ldaupdate/tmp.mat');
end
A1s = cellfun(@(A) reshape(A, prod(p), 1), A1s, 'UniformOutput', 0);
A1s = arrayfun(@(t) cell2mat(A1s(t,:)), 1:trials, 'UniformOutput', 0);
A2s = cellfun(@(A) reshape(A, prod(p), 1), A2s, 'UniformOutput', 0);
A2s = arrayfun(@(t) cell2mat(A2s(t,:)), 1:trials, 'UniformOutput', 0);


disp(' ');