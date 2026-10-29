%PRE_FULL_PSI Precomputation flag
%   If this flag is set, the convolution step (the multiplication with the sparse
%   matrix B) uses (2m+2)^dM precomputed values of the window function, in
%   addition indices of source and target vectors are stored.
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
function f = PRE_FULL_PSI()

f = bitshift(1, 5);
