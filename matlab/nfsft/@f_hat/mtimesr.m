function r = mtimesr(p,q)
% Multiply for f_hat
%
% Copyright 2007 Jens Keiner

if (isa(p,'f_hat') && isnumeric(q))
  p = f_hat(p);
  if (p.N ~= length(q(:))-1)
    error('Dimensions must agree.')
  end
  r = f_hat(p);
  for k = 0:r.N
    r.f_hat(k^2+1:k^2+1+2*k) = q(k+1)*r.f_hat(k^2+1:k^2+1+2*k);
  end
else
  error('Wrong argument types');
end
