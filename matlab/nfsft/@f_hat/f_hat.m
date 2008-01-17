function p = f_hat(a)
% Constructor for f_hat class
%
% Copyright 2007 Jens Keiner

if (nargin == 0)
  p.N = -1;
  p.f_hat = [];
  p = class(p,'f_hat');
elseif (isa(a,'f_hat'))
  p = a;
elseif (isnumeric(a) && isscalar(a) && a >= 0)
  p.N = a;
  p.f_hat = zeros((p.N+1)^2,1);
  p = class(p,'f_hat');
elseif (isnumeric(a) && ndims(a) == 2 && (size(a,1)-1)/2+1 == size(a,2))
  p.N = (size(a,1)-1)/2;
  p.f_hat = zeros((p.N+1)^2,1);
  o = 1;
  for k = 0:p.N
    p.f_hat(o:o+2*k) = a(k*(2*p.N+1)+p.N+1-k:k*(2*p.N+1)+p.N+1+k);
    o = o + 2*k+1;
  end
  p = class(p,'f_hat');
elseif (isnumeric(a) && ndims(a) == 2 && length(a) == 0 || ...
        round(sqrt(length(a(:)))-1) == sqrt(length(a(:)))-1)
  p.N = sqrt(length(a(:)))-1;
  p.f_hat = a(:);
  p = class(p,'f_hat');
else
  error('Invalid argument for f_hat constructor.')
end