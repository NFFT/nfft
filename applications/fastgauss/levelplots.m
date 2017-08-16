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
clear;
n=128;
p=1;
a = 1:1500;
b = -1500:1500;
[A,B] = meshgrid(a,b);
S = A.^2 + B.^2;
Z = exp(-pi^2*A*n^2./(4*p^2*S));
ZZ = 2*exp(-A*(2*p-1)^2/4).*(1+1./p*(2*p-1)*A) + sqrt(pi./sqrt(S))/p.*exp(-pi^2*A*n^2./(4*p^2*S)).*(1 + 2*S*p^2./(n*pi^2*A));

p =-20:2:2;
v = 10.^p;
figure(1);
[C,h]=contour(A,B,ZZ,v);
clabel(C,h);
figure(2);
[C,h]=contour(A,B,Z,v);
clabel(C,h);
