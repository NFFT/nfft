% The following function is based on code by Greg von Winckel, 02/25/2004
% See: TODO Add references.
function [x,w] = lgwt(n,a,b)

n = n-1; n1= n + 1; n2 = n + 2;
% n1 uniform nodes in [-1,1]
xu = linspace(-1,1,n1)';
% initial guess for nodes
y=cos((2*(0:n)'+1)*pi/(2*n+2)) + (0.27/n1)*sin(pi*xu*n/n2);
% Gauss-Legendre Vandermonde matrix
L=zeros(n1,n2);
% derivative of that
Lp=zeros(n1,n2);
% We compute the zeros of the n1th Legendre polynomial using the recursion
% relation and the Newton-Raphson method.
y0=2;
% Iterate until new points are uniformly within epsilon of old points.
while (max(abs(y - y0)) > eps)
  L(:,1) = 1; Lp(:,1) = 0; L(:,2) = y; Lp(:,2) = 1;
  for k = 2:n1
    L(:,k+1) = ((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
  end
  Lp = n2 * (L(:,n1) - y.*L(:,n2))./(1-y.^2); y0 = y; y = y0 - L(:,n2)./Lp;
end
% linear map from[-1,1] to [a,b]
x = (a*(1-y)+b*(1+y))/2;
if (nargout == 2)
  % weights
  w = (b-a)./((1-y.^2).*Lp.^2)*(n2/n1)^2;
end
