%NFSFT_PRECOMPUTE Node-independent precomputation (for FPT)
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
function nfsft_precompute(N,kappa,nfsft_flags,fpt_flags)
if ~exist('kappa','var')
  kappa = 1000;
end

if ~exist('nfsft_flags','var')
  nfsft_flags = 0;
end

if ~exist('fpt_flags','var')
  fpt_flags = 0;
end

nfsftmex('precompute',N,kappa,nfsft_flags,fpt_flags)
