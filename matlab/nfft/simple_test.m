% $Id: simple_test.m 2594 2008-11-14 09:42:46Z keiner $
% 
% Copyright (c) 2005, 2008 Jens Keiner, Stefan Kunis, Daniel Potts
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
% 
%  This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
%  details.
%  
%  You should have received a copy of the GNU General Public License along with
%  this program; if not, write to the Free Software Foundation, Inc., 51
%  Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

% maximum degree (bandwidth)
N = 10;

% number of nodes
M = 9;

% nodes
x=rand(1,M)-0.5;

% Create plan.
plan = nfft_init_1d(N,M);

% Set nodes.
nfft_set_x(plan,x);

% node-dependent precomputation
nfft_precompute_psi(plan);

% Fourier coefficients
f_hat = rand(N,1)+i*rand(N,1);

% Set Fourier coefficients.
nfft_set_f_hat(plan,double(f_hat));

% transform
nfft_trafo(plan);

% function values
f = nfft_get_f(plan)

% finalize plan
nfft_finalize(plan);

A=exp(-2*pi*i*x'*(-N/2:N/2-1));
f2=A*f_hat
