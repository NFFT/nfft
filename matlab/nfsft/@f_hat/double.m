function d = double(p)
% Double conversion for f_hat class
%
% Copyright 2007 Jens Keiner

d = zeros(2*p.N+1,p.N+1);
o = 1;
for k = 0:p.N
  d(k*(2*p.N+1)+p.N+1-k:k*(2*p.N+1)+p.N+1+k) = p.f_hat(o:o+2*k);
  o = o + 2*k+1;
end
