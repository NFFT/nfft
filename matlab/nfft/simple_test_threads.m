%SIMPLE_TEST_THREADS Example program: Basic usage principles for the
%set_num_threads function
%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

nthreads = nfft_get_num_threads();

% degree (bandwidth)
N = 32;

% number of nodes
M = 10000;

% nodes
x=rand(3,M)-0.5;

% Fourier coefficients
f_hat = rand(N^3,1)+1i*rand(N^3,1);

for threads = 1:nthreads
  nfft_set_num_threads(threads);

  % Create plan.
  plan = nfft_init_3d(N,N,N,M);

  % Set nodes.
  nfft_set_x(plan,x);

  % node-dependent precomputation
  nfft_precompute_psi(plan);

  % Set Fourier coefficients.
  nfft_set_f_hat(plan,double(f_hat));

  % transform
  tic
  nfft_trafo(plan);
  t = toc;

  % function values
  f = nfft_get_f(plan);

  % finalize plan
  nfft_finalize(plan);

  fprintf('Threads:%2d,  Time: %1.2e\n',threads,t);
end
