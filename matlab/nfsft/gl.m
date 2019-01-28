%GL Gauss-Legendre interpolatory quadrature rule
%   X = GL(N) generates the (2N+2)*(N+1) Gauss-Legendre nodes and returns a
%   2x[(2N+2)*(N+1)] matrix X containing their spherical coordinates. The first
%   row contains the longitudes in [0,2pi] and the second row the colatitudes in
%   [0,pi].
%
%   [X,W] = GL(N) in addition generates the quadrature weights W. The resulting
%   quadrature rule is exact up to polynomial degree 2*N.
%
%   Example:
%   [X,W] = GL(1) gives
%
%   X =
%        0    1.5708    3.1416    4.7124         0    1.5708    3.1416    4.7124
%   0.9553    0.9553    0.9553    0.9553    2.1863    2.1863    2.1863    2.1863
%
%   W =
%   1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708
%
%   See also cc, healpix, equidist
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
function [x,w] = gl(n)

% correctness conditions
if (nargin ~= 1)
  error('gl: Exactly one input argument required.');
end
if (~isscalar(n) || n < 0)
  error('gl: Input argument must be a non-negative number.');
end
if (nargout > 2)
  error('gl: No more than two output arguments allowed.');
end

% Gauss-Legendre nodes and weights (maybe) on [-1,1]
if (nargout == 2)
  [theta,w] = lgwt(n+1,-1,1);
else
  theta = lgwt(n+1,-1,1);
end
% coordinate transformation to [0,pi] for colatitude
theta = acos(theta);
% equispaced nodes for longitude
phi = 2*pi*(0:2*n+1)/(2*n+2);
% tensor product
[x,y] = meshgrid(theta,phi);
% linearisation to 2x((2*n+2)(n+1)) matrix
x = [y(:),x(:)]';
if (nargout == 2)
  % weights
  w = repmat((pi/(n+1))*w,1,2*n+2)';
  w = w(:)';
end
