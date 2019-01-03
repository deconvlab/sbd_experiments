function x0 = sqclumpx(sz, params)
%SINCLUMPX  Clumped (Bernoulli) sparse activation with square wave sparsity rate.
%   sz (scalar):  sz of activation map.
%   params  (4 nnfloat):  [p0 p1 r1 fq], square wave parameters.

if nargin<2 || isempty(params)
  params = [0.1 0.5 0.2 5];
end

if numel(sz) == 1
  sz = [sz 1];
else
  sz = sz(:)';
end

p0 = params(1);  % firing probability when square wave is off
p1 = params(2);  % firing probability when square wave is on
r1 = params(3);  % on ratio per wave period; if gt 0, wave is on then off, vice versa
fq = params(4);  % number of square wave periods in x0

wvsz = ceil(sz(1)/fq);
sqwv = zeros(wvsz,1);

sqwv(1:round(abs(r1)*wvsz)) = 1;
if r1<0;  sqwv = flipud(sqwv);  end

sqwv = repmat(sqwv, [fq sz(2:end)]);
sqwv = p0 + (p1-p0)*sqwv(1:sz(1),:);

x0 = double(rand(sz) <= sqwv);
end

