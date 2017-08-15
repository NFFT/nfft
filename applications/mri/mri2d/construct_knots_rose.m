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
function [M] = construct_knots_rose( N )

N=ceil(1.5*N);
B=1;
M=N^2;
file = zeros(M,2);

A=0.5;

w=N/64*50;

for b=0:B-1,
for i=1:M/B,
    t=i/(M/B);
    file(b*M/B+i,1) = A*cos(2*pi*w*t)*cos(2*pi*t+b*2*pi/B);

    file(b*M/B+i,2) = A*cos(2*pi*w*t)*sin(2*pi*t+b*2*pi/B);
end
end

% feel free to plot the knots by uncommenting
% plot(file(:,1),file(:,2),'.-');

save knots.dat -ascii file
