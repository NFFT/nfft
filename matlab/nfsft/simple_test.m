%SIMPLE_TEST Example program: Basic usage principles
%   Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts

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
% $Id$

disp(sprintf('Number of threads: %d\n', nfsft_get_num_threads()));

% maximum degree (bandwidth) of spherical harmonics expansions
N = 2;

% threshold (affects accuracy, 1000 is the default)
kappa = 1000;

% precomputation
nfsft_precompute(N,kappa);

% Gauﬂ-Legendre interpolatory quadrature nodes for N. See gl.m
[X,W] = gl(N);

% number of nodes
M = size(X,2);

% Create plan.
plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);

% Set nodes.
nfsft_set_x(plan,X);

% node-dependent precomputation
nfsft_precompute_x(plan);

% Fourier coefficients
fh = f_hat(rand((N+1)*(N+1),1));
display(fh);

% Set Fourier coefficients.
nfsft_set_f_hat(plan,double(fh));

% transform
nfsft_trafo(plan);

% function values
f = nfsft_get_f(plan);
display(f);

% adjoint transform, using quadrature weights tor ecover Fourier coefficients
nfsft_set_f(plan,f.*W');
nfsft_adjoint(plan);

fh2 = f_hat(nfsft_get_f_hat(plan));
display(fh2);

nfsft_finalize(plan);
