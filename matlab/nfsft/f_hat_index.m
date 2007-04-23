function ind = f_hat_index(f_hat,j,k)
% Return indices of Fourier coefficients
%
% Copyright 2007 Jens Keiner

[m,n] = size(f_hat);
N = (m-1)/2;

if (mod(m,2) ~= 1 || N+1 ~= n)
  error('This is not a valid f_hat variable')
end

if (exist('k','var') == 0 || exist('n','var') == 0)
  ind = zeros(1,(N+1)*(N+1));
  j = 1;
  for k = 0:N
    ind(j:j+2*k) = k*(2*N+1) + N + 1 + (-k:k);
    j = j + 2*k+1;
  end
else
  ind = (2*N+1)*j + k + N + 1;
end