%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
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
%   kernel = 'multiquadric', e.g.
%   c kernel parameter
%   m cut-off parameter for NFFT
%   n expansion degree
%   p smoothness
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
system(sprintf('./fastsum_matlab %d %d %d %d %d %d %s %e %e %e',d,N,M,n,m,p,kernel,c,eps_I,eps_B));

%read result from file
f2=load('f.dat');
f=f2(:,1)+i*f2(:,2);

f2=load('f_direct.dat');
f_direct=f2(:,1)+i*f2(:,2);
