% Copyright (c) 2005 The NFFT Authors
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
in=load('output_phantom_nfft_600.dat');

out = [in(:,1) in(:,2)];

tmp=in(:,3)+i*in(:,4);
tmp=tmp.*exp(2*pi*i* (-21*in(:,2) + in(:,1)));
out = [out real(tmp) imag(tmp)];

save output_phantom_nfft_600.dat -ascii out
