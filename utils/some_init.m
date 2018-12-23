function ainit = some_init(solver, p0, a0, x0, calc_grad)
  if nargin < 5;  calc_grad = false;  end
  m = numel(solver.y);

  ainit = solver.y(mod(randi(m) + (1:p0), m) + 1);
  ainit = [zeros(p0-1,1); ainit(:); zeros(p0-1,1)];
  if calc_grad
    ainit = -solver.calc_grad(ainit/norm(ainit(:)));
  end
  ainit = ainit(:)/norm(ainit(:));
end