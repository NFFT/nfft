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
function [M] = construct_knots_radial( N )

A=N*2;
P=410/512*A;
M=A*P;
file=zeros(P*A,2);
for i=1:P
  for j=1:A,
   if mod(i,2) == 0,
     r=(j-1)/A - 1/2;
   else
     r=-(j-1)/A + 1/2;
   end
   file((i-1)*A+j,1)=r*sin((i-1)*pi/P);
   file((i-1)*A+j,2)=r*cos((i-1)*pi/P);
  end
end

% feel free to plot the knots by uncommenting
% plot(file(:,1),file(:,2),'.-');

save knots.dat -ascii file
