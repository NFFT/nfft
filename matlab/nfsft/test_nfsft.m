%TEST_NFSFT Example program: Basic usage principles
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

% maximum degree (bandwidth) of spherical harmonics expansions
N = 512;

% Gauss-Legendre interpolatory quadrature nodes for N. See gl.m
[X,W] = gl(N);

% number of nodes
M = size(X,2);

% Create plan.
p = nfsft(N,M,NFSFT_NORMALIZED);

%% Set nodes.
p.x=X;

% Fourier coefficients
fh = f_hat(rand((N+1)*(N+1),1));
for i=1:N/2
    fh(2*i-1,-2*i+1:2*i-1)=fh(2*i-1,-2*i+1:2*i-1)/2;
end

% Set Fourier coefficients.
p.fhat=fh;

% show Fourier coefficient of degree 2 and order -1
fprintf('fh(2,-1) = %g\n',fh(2,-1))

%% transform
nfsft_trafo(p);

% function values
f = p.f;

%% adjoint transform, using quadrature weights to recover Fourier coefficients
p.f = f.*W';
nfsft_adjoint(p);

fh2 = p.fhat;
fprintf('Error of computation: %g\n',norm(fh2-fh));
