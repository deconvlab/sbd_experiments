function out = soft(in, lam)
    out = sign(in) .* max(abs(in) - lam, 0);
end