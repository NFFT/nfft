%TEST_NFSFT Example program: Advanced usage principles
%Use simple_test.m for the basic usage exampe.

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
fprintf('Number of threads: %d\n', nfsft_get_num_threads());

% maximum degree (bandwidth) of spherical harmonics expansions
N = 128;
fprintf('Performing a fast and direct NFSFT with bandwidth N = %d\n', N);

% threshold (affects accuracy, 1000 is the default)
kappa = 1000;

% random Fourier coefficients
fh = f_hat(rand((N+1)*(N+1),1));

% precomputation
tic
nfsft_precompute(N,kappa);
fprintf('Time of precomputation:    %g seconds\n', toc);

% Gauss-Legendre interpolatory quadrature nodes for N. See gl.m
[X,W] = gl(N);

% number of nodes
M = size(X,2);

% Create plan.
plan = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT));

% Set nodes.
nfsft_set_x(plan,X);

% node-dependent precomputation
nfsft_precompute_x(plan);

% set Fourier coefficients
nfsft_set_f_hat(plan,double(fh));

% direct (slow) transform
tic
ndsft_trafo(plan);
fprintf('Time of direct transform:  %g seconds\n', toc);

% get function values
f_direct = nfsft_get_f(plan);

% fast NFSFT transform
tic
nfsft_trafo(plan);
fprintf('Time of fast transform:    %g seconds\n', toc);

% get function values
f_fast = nfsft_get_f(plan);

% adjoint transform, using quadrature weights to recover Fourier coefficients
nfsft_set_f(plan,f_direct.*W')
tic
ndsft_adjoint(plan);
times.adjoint_direct=toc;
fprintf('Time of direct adjoint:    %g seconds\n', toc);

% get function values
fh2_direct = f_hat(nfsft_get_f_hat(plan));

% Fast adjoint transform
nfsft_set_f(plan,f_direct.*W')
tic
nfsft_adjoint(plan);
times.adjoint=toc;
fprintf('Time of fast adjoint:      %g seconds\n', toc);
fh2_fast = f_hat(nfsft_get_f_hat(plan));

% finalize plan
nfsft_finalize(plan);


% display error
fprintf('Relative error of fast transform:         %g\n', norm(f_fast-f_direct)/norm(f_direct));
fprintf('Relative error of fast adjoint transform: %g\n', norm(fh2_fast-fh2_direct)/norm(fh2_direct));

% destroy all precomputations
nfsft_forget()