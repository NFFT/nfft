%NFST_GET_NUM_THREADS Get number of threads (at most) used for computation
%   Copyright (c) 2009, 2012, 2017 The NFFT Authors

% Copyright (c) 2009, 2012, 2017 The NFFT Authors
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
function nthreads = nfst_get_num_threads()

nthreads = nfstmex('get_num_threads');
