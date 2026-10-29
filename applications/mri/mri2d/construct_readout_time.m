% Copyright (c) 2005, 2009 The NFFT Authors
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <https://www.gnu.org/licenses/>.
%
% constructs a three dimensional Shepp Logan phantom
function [M] = construct_readout_time( M, Te, arms, increment )

out=zeros(M,1);

for j=1:M,
  out(j) =   mod(j-1,M/arms)*increment + Te;
end

%0.00402542373

save readout_time.dat -ascii out
