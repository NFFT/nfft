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
function [] = precompute_weights ( file, M )

input=load(file);

kx = input(1:M,1);
ky = input(1:M,2);

kxy=[kx ky];

% compute the voronoi regions
[V,C] = voronoin(kxy,{'QJ'});

% the surface of the knots is written to area
area = [];

% sum of all surfaces
sum_area = 0;

% the maximum distance two nearest neighbour have
% to get the surface we store max_distance^2
max_distance=0;

% compute the surface of the knots
for j= 1:length(kxy)
  x = V(C{j},1);
  y = V(C{j},2);
  lxy = length(x);
  if(lxy==0) % a knot exists more than one time
    A=0;
  else
    A = abs(sum(0.5*(x([2:lxy 1])-x(:)).* ...
        (y([2:lxy 1]) + y(:))));
  end
  area = [area A];
  min_distance = min((2*(x-kxy(j,1))).^2+(2*(y-kxy(j,2))).^2);
  max_distance = max([max_distance min_distance]);
end

% if the surface of a knot is bigger than max_distance^2
% or isnan or isinf, then take max_distance^2
for j=1:length(area),
  if(isnan(area(j)) | isinf(area(j))| area(j)>max_distance),
    area(j)=max_distance;
  end
  sum_area = sum_area + area(j);
end

% norm the weights
area = area / sum_area;

save weights.dat -ascii area

