%F_HAT Constructor for f_hat class
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
function p = f_hat(a)

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
elseif (isnumeric(a) && ndims(a) == 2 && isempty(a) || ...
        round(sqrt(numel(a))-1) == sqrt(numel(a))-1)
  p.N = sqrt(numel(a))-1;
  p.f_hat = a(:);
  p = class(p,'f_hat');
else
  error('Invalid argument for f_hat constructor.')
end
