%CC Clenshaw-Curtis interpolatory quadrature rule
%   X = CC(N) generates the (2N+2)*(2*N-1)+2 Clenshaw-Curtis nodes and returns a
%   2x[(2N+2)*(2*N-1)+2] matrix X containing their spherical coordinates. The
%   first row contains the longitudes in [0,2pi] and the second row the
%   colatitudes in [0,pi].
%
%   [X,W] = CC(N) in addition generates the quadrature weights W. The resulting
%   quadrature rule is exact up to polynomial degree 2*N.
%
%   For the special case N = 0, the number of nodes is 2.
%
%   Nota bene: In latitudinal direction, the nodes are based on Chebyshev zeros
%   in (-1,1). Similar quadrature rules based on Chebyshev extrema in [-1,1]
%   are also often denoted Clenshaw-Curtis quadrature rules. See the references.
%
%   Example:
%   [X,W] = CC(1) gives
%
%   X =
%        0         0    1.5708    3.1416    4.7124         0
%        0    1.5708    1.5708    1.5708    1.5708    3.1416
%
%   W =
%   2.0944    2.0944    2.0944    2.0944    2.0944    2.0944
%
%   See also gl, healpix, equidist
%
%   References
%   TODO Add references.
%
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
function [x,w] = cc(n)

% correctness conditions
if (nargin ~= 1)
  error('cc: Exactly one input argument required.');
end
if (~isscalar(n) || n < 0)
  error('cc: Input argument must be a non-negative number.');
end
if (nargout > 2)
  error('cc: No more than two output arguments allowed.');
end

% Gauss-Legendre nodes and weights (maybe) on [-1,1]
if (nargout == 2)
  [theta,w] = fclencurt(2*n+1,-1,1);
else
  theta = fclencurt(2*n+1,-1,1);
end
% coordinate transformation to [0,pi] for colatitude
theta = acos(theta);
% equispaced nodes for longitude
phi = 2*pi*(0:2*n+1)/(2*n+2);
% tensor product
[x,y] = meshgrid(theta(2:end-1),phi);
% linearisation to 2x((2*n+2)(2*n-1))+2 matrix
x = [y(:),x(:)]';
x = [[0;0],x,[0;pi]];
if (nargout == 2)
  % weights
  w0 = 2*pi*w(1);
  w = repmat((pi/(n+1))*w(2:end-1),1,2*n+2)';
  w = w(:)';
  w = [w0,w,w0];
end

% The following function is based on code by Greg von Winckel, 02/12/2005
% See: ...
function [x,w] = fclencurt(n1,a,b)
n = n1-1; bma = b - a;
c = zeros(n1,2);
c(1:2:n1,1) = (2./[1 1-(2:2:n).^2 ])'; c(2,2)=1;
f = real(ifft([c(1:n1,:);c(n:-1:2,:)]));
w = bma*([f(1,1); 2*f(2:n,1); f(n1,1)])/2;
x = 0.5*((b+a)+n*bma*f(1:n1,2));

