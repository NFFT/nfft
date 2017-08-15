%FASTSUM_SIMPLE Compute fast and, optionally, direct fast summation.
%
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% Computes the sums
% 
%   f(y_j) = sum_{k=1}^N alpha_k kernel(x_k-y_j)   (j=1:M)
% 
% Calling:
%   f = fastsum(x,alpha,y,kernel,c,m,n,p,eps_I,eps_B)
%   [f,f_direct] = fastsum(x,alpha,y,kernel,c,m,n,p,eps_I,eps_B)
% 
% Output:
% f           fast summation
% f_direct    computed with direct summation
% 
% Input:
% size(x)     = [N,d]   source knots in d-ball with radius 1/4-eps_b/2
% size(alpha) = [N,1]   source coefficients (complex)
% size(y)     = [M,d]   target knots in d-ball with radius 1/4-eps_b/2
% size(f)     = [N,1]   target evaluations (complex)
% kernel = 'multiquadric', etc. (see options below)
% c       kernel parameter
% n       expansion degree
% m       cut-off parameter for NFFT
% p       degree of smoothness of regularization
% eps_I   inner boundary
% eps_B   outer boundary
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
function [f,f_direct]=fastsum_simple(x,alpha,y,kernel,c,m,n,p,eps_I,eps_B...
    ,nn_oversampled,flags)

nargoutchk(1, 2)
d = size(x,2);

if(nargin<12)
    flags = 0;
end
if (nargin<11)
    nn_oversampled = 2*n;
end

plan=fastsum_init(d,kernel,c,flags,n,p,eps_I,eps_B);
fastsum_set_x_alpha(plan,x,alpha,nn_oversampled,m)
fastsum_set_y(plan,y,nn_oversampled,m)

if(nargout==2)
    fastsum_trafo_direct(plan)   % direct computation
    f_direct = fastsum_get_f(plan);
end

fastsum_trafo(plan)         % fast computation
f = fastsum_get_f(plan);
fastsum_finalize(plan)
end
