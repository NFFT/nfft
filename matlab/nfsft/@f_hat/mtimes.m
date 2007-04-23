function r = mtimes(p,q)
% Multiply for f_hat
%
% Copyright 2007 Jens Keiner

p = f_hat(p);
q = f_hat(q);
if (p.N ~= q.N)
  error('Dimensions must agree.')
end
r = f_hat(p.f_hat .* q.f_hat);
