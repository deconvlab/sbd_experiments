function ainit = data_ainit(solver, p0, calc_grad)
  assert(~isempty(solver.y), 'y has not been initialized.');
  if nargin < 3;  calc_grad = false;  end
  m = numel(solver.y);

  ainit = solver.y(mod(randi(m) + (1:p0), m) + 1);
  ainit = [zeros(p0-1,1); ainit(:); zeros(p0-1,1)];
  if calc_grad
    ainit = -solver.calc_grad(ainit/norm(ainit(:)));
  end
  ainit = ainit(:)/norm(ainit(:));
end