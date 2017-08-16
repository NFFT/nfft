%SIMPLE_TEST Example program: Basic usage principles
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
fprintf('Number of threads: %d\n', nfct_get_num_threads());
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('A simple one dimensional example');

% maximum degree (bandwidth)
N = 14;

% number of nodes
M = 19;

% nodes
x=0.5*rand(1,M);

% Create plan.
plan = nfct_init_1d(N,M);

% Set nodes.
nfct_set_x(plan,x);

% node-dependent precomputation
%nfct_precompute_psi(plan);

% Fourier coefficients
f_hat = rand(N,1);

% Set Fourier coefficients.
nfct_set_f_hat(plan,double(f_hat));

% transform
nfct_trafo(plan);

% function values
f = nfct_get_f(plan);

nfct_adjoint(plan);

f_hat_adjoint = nfct_get_f_hat(plan);

% finalize plan
nfct_finalize(plan);

A=cos(2*pi*x'*(0:N-1));
f2 = A*f_hat;

f_hat_adjoint2 = A'*f;

error_vector = f-f2;
error_linfl1 = norm(f-f2,inf)/norm(f_hat,1);
fprintf('error trafo: %.1e\n', error_linfl1);

error_vector_adjoint = f_hat_adjoint-f_hat_adjoint2;
error_linfl1_adjoint = norm(f_hat_adjoint-f_hat_adjoint2,inf)/norm(f,1);
fprintf('error adjoint: %.1e\n', error_linfl1_adjoint);
