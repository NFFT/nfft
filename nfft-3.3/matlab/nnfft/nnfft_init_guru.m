%NNFFT_INIT_GURU Initialise plans, no error handling
%   Matlab might run into a segmentation violation for wrong parameters
%
%   nnfft_init_guru(d,N1,...,Nd,M,n1,...,nd,m,flags)
%
%   d            spatial dimension
%   N1,...,Nd    bandwidths
%   M            number of nodes
%   n1,...,nd    fft lengths
%   m            cut-off parameter
%   flags   PRE_PHI_HUT | {FG_PSI, PRE_LIN_PSI, PRE_FG_PSI, PRE_PSI,
%	             PRE_FULL_PSI} | FFT_OUT_OF_PLACE
%
%   Copyright (c) 2002, 2015 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2015 Jens Keiner, Stefan Kunis, Daniel Potts
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
%
% $Id: nnfft_init_guru.m 4086 2014-12-11 19:09:50Z kunis $
function p = nnfft_init_guru(varargin)

p = nnfftmex('init_guru',varargin);
