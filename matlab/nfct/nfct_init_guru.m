%NFCT_INIT_GURU Initialise plans, no error handling
%   Matlab might run into a segmentation violation for wrong parameters
%
%   nfct_init_guru(d,N1,...,Nd,M,n1,...,nd,m,flags,fftw_flags)
%
%   d            spatial dimension
%   N1,...,Nd    bandwidths
%   M            number of nodes
%   n1,...,nd    fft lengths
%   m            cut-off parameter
%   flags   PRE_PHI_HUT | {FG_PSI, PRE_LIN_PSI, PRE_FG_PSI, PRE_PSI,
%	             PRE_FULL_PSI} | FFT_OUT_OF_PLACE
%   fftw_flags   {FFTW_ESTIMATE, FFTW_MEASURE}
%
%   Copyright (c) 2009, 2013, 2017 The NFFT Authors

% Copyright (c) 2009, 2013, 2017 The NFFT Authors
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
function p = nfct_init_guru(varargin)

p = nfctmex('init_guru',varargin);
