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
N = 2000;       % number of source knots
M = 2000;       % number of target knots
kernel = 'multiquadric';
c = 1/sqrt(N);  % kernel parameter
p = 3;          % degree of smoothness of regularization
flags = 0;      % flags (could be EXACT_NEARFIELD or NEARFIELD_BOXES)
n = 156;        % expansion degree
eps_I = p/n;    % inner boundary
eps_B = 1/16;   % outer boundary
m = p;          % cut-off parameter for NFFT
nn_oversampled=2*n; % oversampling factor for NFFT

%% random source nodes in circle of radius 0.25-eps_B/2
r = sqrt(rand(N,1))*(0.25-eps_B/2);
phi = rand(N,1)*2*pi;
x = [r.*cos(phi) r.*sin(phi)];
% random coefficients
alpha = rand(N,1) + 1i*rand(N,1);
% random target nodes in circle of radius 0.25-eps_B/2
r = sqrt(rand(M,1))*(0.25-eps_B/2);
phi = rand(M,1)*2*pi;
y = [r.*cos(phi) r.*sin(phi)];

%% Perform fastsum
plan=fastsum_init(d,kernel,c,flags,n,p,eps_I,eps_B);
fastsum_set_x(plan,x,nn_oversampled,m)
fastsum_set_alpha(plan,alpha)
fastsum_set_y(plan,y,nn_oversampled,m)

fastsum_trafo_direct(plan)   % direct computation
f_dir = fastsum_get_f(plan);

fastsum_trafo(plan)         % fast computation
f = fastsum_get_f(plan);
fastsum_finalize(plan)

%% fastsum class
a=fastsum(d,kernel,c,flags,n,p,eps_I,eps_B,nn_oversampled,m);
a.x=x;
a.alpha=alpha;
a.y=y;
fastsum_trafo(a);
disp(max(a.f-f))

%% plot source and target evaluations
figure(1)
scatter(x(:,1),x(:,2),120,real(alpha),'.')
colorbar
figure(2)
scatter(y(:,1),y(:,2),120,abs(f./f_dir-1),'.')
colorbar
