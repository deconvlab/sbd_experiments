function ainit = some_init(solver, p0, a0, x0)
  m = numel(solver.y);

  ainit = solver.y(mod(randi(m) + (1:p0), m) + 1);
  ainit = [zeros(p0-1,1); ainit(:); zeros(p0-1,1)];
  ainit = -solver.calc_grad(ainit/norm(ainit(:)));
  ainit = ainit(:)/norm(ainit(:));
end