% Modified third-party code, based on lgwt.m by Greg von Winckel, 02/25/2004,
% from MATLAB Central File Exchange, file ID 4540,
% "Legendre-Gauss Quadrature Weights and Nodes",
% https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
% It is not covered by the NFFT license, but by the following license:
%
% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
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
