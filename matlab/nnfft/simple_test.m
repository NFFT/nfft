%SIMPLE_TEST Example program: Basic usage principles
%
%   Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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
% $Id: simple_test.m 3776 2012-06-03 13:29:25Z keiner $

disp(sprintf('Number of threads: %d\n', nnfft_get_num_threads()));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('A simple one dimensional example');

% maximum degree (bandwidth)
N = 10;

% number of nodes
M = 3;

% nodes
x=rand(1,M)-0.5;

% Create plan.
plan = nnfft_init_1d(N,M);
 %plan = nnfft_init_guru(1,N,M,N,2*N,2,bitor(PRE_PHI_HUT,PRE_PSI));

% Set nodes.
nnfft_set_x(plan,x);
nnfft_get_x(plan);

% node-dependent precomputation
nnfft_precompute_psi(plan);

% Fourier coefficients
f_hat = rand(N,1)+i*rand(N,1);

% Set Fourier coefficients.
nnfft_set_f_hat(plan,double(f_hat));
nnfft_get_f_hat(plan);

% transform
nnfft_trafo(plan);

% function values
f = nnfft_get_f(plan)

% finalize plan
nnfft_finalize(plan);

%A=exp(-2*pi*i*x'*(-N/2:N/2-1));
%f2=A*f_hat

%error_vector=f-f2
%error_linfl1=norm(f-f2,inf)/norm(f_hat,1)

