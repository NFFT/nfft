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
function [M] = construct_knots_linogram ( N,Z )

l=((-N/2):(N/2-1));
k=((-N):(N-1))';

x=[reshape(k*l/N^2,2*N^2,1) kron(ones(N,1),k/2/N);
   kron(ones(N,1),k/2/N) reshape(-k*l/N^2,2*N^2,1)];

file=[];
for z=0:Z-1,
  file=[file; x ones(size(x,1),1)*(z/Z-0.5)];
end

M=size(file,1);

% feel free to plot the knots by uncommenting
% plot(file(1:M/Z,1),file(1:M/Z,2),'.');


save knots.dat -ascii file

