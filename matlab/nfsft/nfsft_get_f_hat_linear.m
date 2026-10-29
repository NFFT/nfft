%NFSFT_GET_F_HAT_LINEAR Get Fourier coefficients in a linear vector from plan
%   Copyright (c) 2008-2009 The NFFT Authors

% Copyright (c) 2008-2009 The NFFT Authors
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <https://www.gnu.org/licenses/>.
function f_hat = nfsft_get_f_hat_linear(p)
f_hat = nfsftmex('get_f_hat_linear',p);
