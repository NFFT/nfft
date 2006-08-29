function writeWeights(file,m)
% $Id$
%
% generateNodes - Generate Gauss-Legendre quadrature weights
%
% Copyright (C) 2006 Jens Keiner
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

for k = m
  % Generate Gauss-Legendre nodes in co-latitudinal direction.
  [theta,w] = lgwt(k+1,-1,1);
  theta = (1/(2*pi))*acos(theta);

  % Write data to file.
  fprintf(file,'\n');
  fprintf(file,'\n');
  fprintf(file,'%.30f\n',w);
  fprintf(file,'\n');
  fprintf(file,'%.30f\n',theta);
end
