%MTIMESR Multiply for f_hat
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
function r = mtimesr(p,q)
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
