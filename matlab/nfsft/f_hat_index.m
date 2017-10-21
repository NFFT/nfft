%F_HAT_INDEX Return indices of Fourier coefficients
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
function ind = f_hat_index(f_hat,j,k)
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
