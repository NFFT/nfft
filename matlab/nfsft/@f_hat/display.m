function display(p)
% Display function for f_hat class
%
% Copyright 2007 Jens Keiner

if (p.N == -1)
  disp([inputname(1) ' = empty']);
else
  f_hat_m = zeros(2*p.N+1,p.N+1);
  o = 1;
  for k = 0:p.N
    f_hat_m(k*(2*p.N+1)+p.N+1-k:k*(2*p.N+1)+p.N+1+k) = p.f_hat(o:o+2*k);
    o = o + 2*k+1;
  end
  disp([inputname(1) ' =']);
  disp(f_hat_m);
end