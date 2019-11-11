%SIMPLE_TEST Example program: Basic usage principles
% Computes the sums
% f(theta_d,phi_d) = sum_{k=0}^N sum_{n=-k}^k f_hat(k,n) Y_k^n(theta_d,phi_d)
% in terms of spherical harmonics Y_k^n of degree k on a set of arbitrary
% nodes (theta_d,phi_d) for d=1..M, in spherical coordinates 
% 0 <= phi_d <2*pi and 0 <= theta_d <= pi.
% We also recover the original data f from the Fourier coefficients f_hat
% via an exact quadrature rule.

% Copyright (c) 2002, 2019 Jens Keiner, Stefan Kunis, Daniel Potts
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
N = 64;

% Gauss-Legendre exact quadrature nodes for degree N. See gl.m
[X,W] = gl(N);
% Clenshaw-Curtis exact quadrature rule is available in cc.m
% Further exact quadrature rules are available at
% https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/
% They can be imported with either the command
%    X = importdata('Design_1002000_1000_random.dat')';
% or
%    X = importdata('N124_M7812_Ico.dat',' ',2);
%    X = [atan2(X.data(:,1),X.data(:,2))';acos(X.data(:,3))'];
% and then
%    W = ones(1,length(X))/length(X)*4*pi;
% These quadrature rules have to be at least twice the size of N,
% e.g. for 'N124_M7812_Ico.dat' this means N<=64.

% number of nodes
M = size(X,2);

% Create plan of class NFSFT.
plan = nfsft(N,M,NFSFT_NORMALIZED);

% Set nodes in spherical coordinates x = [phi; theta]
plan.x = X; 

% random Fourier coefficients
fh = f_hat(rand((N+1)*(N+1),1));

% Set Fourier coefficients.
plan.fhat = fh;

% show Fourier coefficient of degree 2 and order -1
fprintf('fh(2,-1) = %g\n',fh(2,-1))

% NFSFT transform
nfsft_trafo(plan);

% function values
f = plan.f;

% adjoint transform, using quadrature weights to recover Fourier coefficients
% the adjoint is only the inverse with the right quadrature weights W
plan.f = f.*W';
nfsft_adjoint(plan);

fh2 = plan.fhat;
fprintf('Relative error of reconstructed f_hat: %g\n',norm(fh2-fh)/norm(fh));
