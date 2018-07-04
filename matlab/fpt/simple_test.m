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

% The computations performed by the fpt are based on Potts, D. Fast algorithms
% for discrete polynomial transforms on arbitrary grids. Linear Algebra Appl.
% 366, 353-370., 2003  (https://www-user.tu-chemnitz.de/~potts/paper/ndct.pdf)

% This example shows the use of the fast polynomial transform to compute
% the Chebyshev coefficients from a polynomial given by Legendre coefficients,
%
%   f(x) = a_0 P_0(x) + a_1 P_1(x) + ... + a_k_end P_k_end(x).                   (1)
t = 8;
N_max = 2^t; % maximal polynomial degree
k_end = N_max;
k_start = 0;

fprintf('Initializing fast polynomial transform for degree up to N = %d\n', N_max);
% An fpt_set is a data structure that contains precomputed data for a polynomial
% transform. The first parameter (t) is the exponent of the maximum transform size
% desired (2^t), i.e., t = 3 means that N_max in (1) can be at most N = 8.
% The second parameter are flags.
fpt_set = fpt_init(t,0);

% Three-term recurrence coefficients for Legendre polynomials
% alpha(1) and beta(1) are not referenced.
% gamma(1) contains the value of P_0(x) (which is a constant).
% Note that the length of these vectors has to be N_max+2
alpha = [0 ((2*(0:N_max)+1)./((0:N_max)+1))]';
beta = zeros(N_max+2,1);
gamma = [1 -((0:N_max)./((0:N_max)+1))]';

fprintf('Precomputations for Legendre polynomials\n');
% The function fpt_precompute actually does the precomputation for a transform.
% It needs arrays alpha, beta, and gamma, containing the three-term recurrence
% coefficients, here of the Legendre polynomials. The format is explained above.
% The fifth parameter (k_start) is where the index in the linear combination (1)
% starts, here k_start=0.
fpt_precompute(fpt_set,alpha,beta,gamma,k_start);

% random Fourier coefficients of the length at most N_max
a = 2*(rand(k_end+1-k_start,1) + 1i*rand(k_end+1-k_start,1)) - 1;
fprintf('l1-norm of coefficients of Legendre polynomial: %.3e\n\n', norm(a,1));

fprintf('Converting coefficients of Legendre polynomial of degree N = %d into Chebyshev basis\n\n', k_end);
% fast polynomial transform for the polynomial degree k_end and flags 0.
b = fpt_trafo(fpt_set,a,k_end,0);
fprintf('l1-norm of coefficients of Chebyshev polynomial: %.3e\n\n', norm(b,1));

% Comparison with dpt (slow direct version of fpt)

b_direct = fpt_trafo_direct(fpt_set,a,k_end,0);
fprintf('Comparison of Chebyshev coefficients fpt vs. fpt_direct (dpt),\nmax abs error = %.3e\n\n', norm(b-b_direct,'inf'));

% cleanup
fpt_finalize(fpt_set);

% Comparison with Legendre polynomial

nodes = 2 * rand(10*(k_end+1),1) - 1;

% first evaluate Legendre polynomial at random nodes using Clenshaw algorithm:
za = [zeros(k_start,1); a];
fct_legendre_clenshaw = eval_clenshaw(alpha, beta, gamma, za, nodes);

% then evaluate Chebyshev polynomial at random nodes using Clenshaw algorithm:
alpha_chebyshev = [0 1 repmat(2,1,N_max)]';
beta_chebyshev = zeros(N_max+2,1);
gamma_chebyshev = [1 repmat(-1,1,N_max)]';
fct_chebyshev = eval_clenshaw(alpha_chebyshev, beta_chebyshev, gamma_chebyshev, b, nodes);
% for fast version, please see non-equispaced fast cosine transform (nfct).

fprintf('Comparison of input Legendre polynomial\n');
fprintf('and output Chebyshev polynomial evaluated at %d random nodes:\n', length(nodes));
fprintf('Chebyshev poly.(Clenshaw) vs. Legendre(Clenshaw),     max abs error = %.3e\n', norm(fct_chebyshev-fct_legendre_clenshaw,'inf'));

% evaluate Chebyshev polynomial at random nodes using direct evaluation:
fct_chebyshev_direct = zeros(length(nodes),1);
for k=0:k_end
  fct_chebyshev_direct = fct_chebyshev_direct + b(k+1) * cos(k * acos(nodes));
end
fprintf('Chebyshev poly.(direct eval.) vs. Legendre(Clenshaw), max abs error = %.3e\n', norm(fct_chebyshev_direct-fct_legendre_clenshaw,'inf'));

% Comparison with Matlab legendre function
f_direct = zeros(size(nodes));
for k=0:k_end
  L = legendre(k, nodes);
  f_direct = f_direct + za(k+1)*L(1,:).';
end
fprintf('Chebyshev poly.(Clenshaw) vs. Legendre(direct eval.), max abs error = %.3e\n', norm(fct_chebyshev-f_direct,'inf'));
