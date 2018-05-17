dot1 = NaN(trials, maxit);  dot2 = dot1;
for t = 1:trials
    dot1(t,:) = A1s{t}(:,end)' * A1s{t};
    dot2(t,:) = A1s{t}(:,end)' * A2s{t};
end