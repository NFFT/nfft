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
function [M] = construct_knots_spiral ( N,arms )
M=N^2;
file = zeros(M,2);

A=0.5;

w=N/64*50;

hold on;

for b=0:arms-1,
  for i=1:M/arms,
    t=((i-1)/(M/arms))^(1/2);
    file(b*M/arms+i,1) = A*t*cos(2*pi*(w*t+b/arms));

    file(b*M/arms+i,2) = A*t*sin(2*pi*(w*t+b/arms));
  end
  % plot(file(b*M/arms+1:(b+1)*M/arms,1),file(b*M/arms+1:(b+1)*M/arms,2),'-');
end
hold off;

% feel free to plot the knots by uncommenting
% plot(file(:,1),file(:,2),'.');

save knots.dat -ascii file
