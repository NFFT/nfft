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
% Computes the sums
% 
%   f(y_j) = sum_{k=1}^N alpha_k kernel(x_k-y_j)   (j=1:M).
% 
% size(x)     = [N,d]   source knots in d-ball with radius 1/4-eps_b/2
% size(alpha) = [N,1]   source coefficients (complex)
% size(y)     = [M,d]   target knots in d-ball with radius 1/4-eps_b/2
% size(f)     = [N,1]   target evaluations (complex)
% 
% Kernel functions:
% 'gaussian'                K(x) = EXP(-x^2/c^2) 
% 'multiquadric'            K(x) = SQRT(x^2+c^2)
% 'inverse_multiquadric'    K(x) = 1/SQRT(x^2+c^2)
% 'logarithm'               K(x) = LOG |x|
% 'thinplate_spline'        K(x) = x^2 LOG |x|
% 'one_over_square'         K(x) = 1/x^2
% 'one_over_modulus'        K(x) = 1/|x|
% 'one_over_x'              K(x) = 1/x
% 'inverse_multiquadric3'   K(x) = 1/SQRT(x^2+c^2)^3
% 'sinc_kernel'             K(x) = SIN(cx)/x
% 'cosc'                    K(x) = COS(cx)/x
% 'cot'                     K(x) = cot(cx)
% 'one_over_cube'           K(x) = 1/x^3
% 'log_sin'                 K(x) = LOG(|SIN(cx)|)
% 'laplacian_rbf'           K(x)=EXP(-|x|/c)

%% Initialize parameters
d = 2;          % number of dimensions
N = 20000;      % number of source knots
M = 20000;      % number of target knots
kernel = 'multiquadric';
c = 1/sqrt(N);  % kernel parameter
p = 6;          % degree of smoothness of regularization
flags = 0;      % flags (could be EXACT_NEARFIELD or NEARFIELD_BOXES)
n = 256;        % bandwidth in frequency domain for NFFT
eps_I = p/n;    % inner boundary, nearfield radius
eps_B = max(1/16,p/n); % outer boundary
m = p;          % window cut-off parameter for NFFT
nn_oversampled = 2*n; % oversampled bandwidth in frequency domain for NFFT
%% random source nodes in circle of radius 0.25-eps_B/2
x = [];
while size(x,1) < N
  x_new = rand(N,d) * (0.5-eps_B) - (0.25-eps_B/2);
  x = [x; x_new(sum(x_new.^2,2) < (0.25-eps_B/2)^2, :)];
  clear x_new;
end
x = x(1:N,:);
% random coefficients from [0,1] + [0,1] i
alpha = rand(N,1) + 1i*rand(N,1);
% random target nodes in circle of radius 0.25-eps_B/2
y = [];
while size(y,1) < N
  y_new = rand(M,d) * (0.5-eps_B) - (0.25-eps_B/2);
  y = [y; y_new(sum(y_new.^2,2) < (0.25-eps_B/2)^2, :)];
  clear y_new;
end
y = y(1:M, :);

%% Perform fastsum via class interface
fprintf('fastsum d=%d, N=%d source nodes, M=%d target nodes,\n', d, N, M);
fprintf('        kernel=%s, parameter c=%.3g, fastsum flags=%d,\n', kernel, c, flags);
fprintf('        NFFT bandwidth n=%d, regularization smoothness p=%d,\n', n, p);
fprintf('        inner regularization eps_I=%.6g,\n', eps_I);
fprintf('        outer regularization eps_B=%.6g,\n', eps_B);
fprintf('        nn_oversampled=%d, NFFT window cut-off m=%d\n', nn_oversampled, m);
fprintf('number of threads: %d\n', fastsum_get_num_threads);

tic
fastsum_instance = fastsum(d,kernel,c,flags,n,p,eps_I,eps_B,nn_oversampled,m);
time_init = toc;

tic
fastsum_instance.x = x;
fastsum_instance.alpha = alpha;
fastsum_instance.y = y;
time_precompute = toc;

tic
fastsum_trafo(fastsum_instance);
f_fast = fastsum_instance.f;
time_trafo = toc;

%% Compute reference results via direct computation
tic
fastsum_trafo_direct(fastsum_instance);
f_direct = fastsum_instance.f;
time_direct = toc;

%% Compare results
fprintf('time direct computation: %.3e\n', time_direct);
fprintf('time fastsum init:         %.3e\n', time_init);
fprintf('time fastsum precompute:   %.3e\n', time_precompute);
fprintf('time fastsum computations: %.3e\n', time_trafo);
fprintf('maximal absolute error: %.3e\n', max(abs(f_fast-f_direct)));
fprintf('maximal relative error: %.3e\n', max(abs(f_fast-f_direct)./abs(f_direct)));
fprintf('relative l_2 error:     %.3e\n', norm(f_fast-f_direct)/norm(f_direct));
