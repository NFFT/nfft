%PRE_LIN_PSI Precomputation flag
%   If this flag is set, the convolution step (the multiplication with the sparse
%   matrix B) uses linear interpolation from a lookup table of equispaced samples
%   of the window function instead of exact values of the window function.
%
%   Copyright (c) 2009 The NFFT Authors

% Copyright (c) 2009 The NFFT Authors
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
function f = PRE_LIN_PSI()

f = bitshift(1, 2);
