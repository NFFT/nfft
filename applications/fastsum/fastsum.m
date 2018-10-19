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
function [f,f_direct]=fastsum(x,alpha,y,kernel,c,m,n,p,eps_I,eps_B)

% f=fastsum(x,alpha,y,kernel,c,m,n,p)
%
%   Computes the sums
%
%     f(y_j) = sum_{k=1}^N alpha_k kernel(x_k-y_j)   (j=1:M)
%
%   by calling C-program with the fast NFFT-based algorithm.
%
%   size(f) = [N,1] (complex)
%   size(x) = [N,d]
%   size(alpha) = [N,1] (complex)
%   size(y)=[M,d]
%   kernel = 'multiquadric', etc. (see options below)
%   c kernel parameter
%   m cut-off parameter for NFFT
%   n expansion degree
%   p smoothness
% 
%   Kernel functions:
%   'gaussian'                K(x) = EXP(-x^2/c^2) 
%   'multiquadric'            K(x) = SQRT(x^2+c^2)
%   'inverse_multiquadric'    K(x) = 1/SQRT(x^2+c^2)
%   'logarithm'               K(x) = LOG |x|
%   'thinplate_spline'        K(x) = x^2 LOG |x|
%   'one_over_square'         K(x) = 1/x^2
%   'one_over_modulus'        K(x) = 1/|x|
%   'one_over_x'              K(x) = 1/x
%   'inverse_multiquadric3'   K(x) = 1/SQRT(x^2+c^2)^3
%   'sinc_kernel'             K(x) = SIN(cx)/x
%   'cosc'                    K(x) = COS(cx)/x
%   'cot'                     K(x) = cot(cx)
%   'one_over_cube'           K(x) = 1/x^3
%   'laplacian_rbf'           K(x)=EXP(-|x|/c)
%
% Markus Fenn, 2006.

[N,d]=size(x);
[M,d]=size(y);

%write input to file
save -ascii -double x.dat x
alpha2=[real(alpha) imag(alpha)];
save -ascii -double alpha.dat alpha2
save -ascii -double y.dat y

%execute C-program for fast summation
if ispc
    cmd='fastsum_matlab.exe';
else 
    cmd='./fastsum_matlab';
end
system(sprintf('%s %d %d %d %d %d %d %s %.16g %.16g %.16g',cmd,d,N,M,n,m,p,kernel,c,eps_I,eps_B));

%read result from file
f2=load('f.dat');
f=f2(:,1)+i*f2(:,2);

f2=load('f_direct.dat');
f_direct=f2(:,1)+i*f2(:,2);
