%SIMPLE_TEST Example program: Basic usage principles
%
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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


% This example shows the use of the fast polynomial transform to compute
% the Chebyshev coefficients from a polynomial given by Legendre coefficients,
%
%   f(x) = a_0 P_0(x) + a_1 P_1(x) + ... + a_N P_N(x).                   (1)
t = 3;
N = 2^t; % maximal polynomial degree

% An fpt_set is a data structure that contains precomputed data for a number
% of different polynomial transforms. The first parameter (t) is the exponent
% of the maximum transform size desired (2^t), i.e., t = 3 means that N in (1)
% can be at most N = 8. The third parameter are flags.
fpt_set = fpt_init(t,0);

% Three-term recurrence coefficients for Legendre polynomials
% alpha(1) and beta(1) are not referenced.
% gamma(1) contains the value of P_0(x) (which is a constant).
% Note that the length of these vectors has to be \le N+2
alpha = [0 ((2*(0:N)+1)./((0:N)+1))]';
beta = zeros(N+2,1);
gamma = [1 -((0:N)./((0:N)+1))]';

% The function fpt_repcompute actually does the precomputation for a single
% transform. It needs arrays alpha, beta, and gamma, containing the three-
% term recurrence coefficients, here of the Legendre polynomials. The format
% is explained above. The fith parameter (k_start) is where the index in the
% linear combination (1) starts, here k_start=0.
fpt_precompute(fpt_set,alpha,beta,gamma,0);

% random Fourier coefficients of the same length as alpha, beta and gamma - 1
a = 2*rand(N+1,1)-1;

% fast polynomial transform for the polynomial degree N and flags 0.
b = fpt_trafo(fpt_set,a,N,0);

% cleanup
fpt_finalize(fpt_set);
