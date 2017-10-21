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
% Simple test program for fast NFFT-based summation algorithm.
% Markus Fenn, 2006.

%  d=2;
%  N=4000;
%  M=4000;
%  n=128;
%  m=4;
%  p=3;
%  kernel='multiquadric';
%  c=1/sqrt(N);
%
%  system(sprintf('./fastsum_test %d %d %d %d %d %d %s %e',d,N,M,n,m,p,kernel,c));

N = 2000;
M = 2000;
kernel = 'multiquadric';
c = 1/sqrt(N);
m = 4;
p = 3;
n = 156;
eps_I = p/n;
eps_B = 1/16;

%random source nodes in circle of radius 0.25-eps_B/2
r = sqrt(rand(N,1))*(0.25-eps_B/2);
phi = rand(N,1)*2*pi;
x = [r.*cos(phi) r.*sin(phi)];

%random coefficients
alpha = rand(N,1)+i*rand(N,1);

%random target nodes in circle of radius 0.25-eps_B/2
r = sqrt(rand(M,1))*(0.25-eps_B/2);
phi = rand(M,1)*2*pi;
y = [r.*cos(phi) r.*sin(phi)];

%fast NFFT-based summation
[f,f_direct] = fastsum(x,alpha,y,kernel,c,m,n,p,eps_I,eps_B);
