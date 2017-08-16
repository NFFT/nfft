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
function [M] = construct_knots_linogram ( N )

l=((-N/2):(N/2-1));
k=((-N):(N-1))';

x=[reshape(k*l/N^2,2*N^2,1) kron(ones(N,1),k/2/N);
   kron(ones(N,1),k/2/N) reshape(-k*l/N^2,2*N^2,1)];

% feel free to plot the knots by uncommenting
% plot(x(:,1),x(:,2),'r.')

save knots.dat -ascii x

M=size(x,1);
