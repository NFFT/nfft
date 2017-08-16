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
function [] = visualize_data ( file, N, num_fig, caption )

load output_real.dat
load output_imag.dat

A=reshape(abs(output_real+i*output_imag),N,N);

for k=1:N,
  for l=1:N,
    if sqrt((l-N/2)^2+(k-N/2)^2)>N/2,
      A(k,l)=0;
    end
  end
end

% plot the two dimensional phantom
figure(2*num_fig-1)
imagesc(A, [0 1.2]);
colormap(gray(256));
colorbar;
title(caption);
file_out =[file '.jpg'];
print('-djpeg',file_out);

% plot the N/2 row
figure(2*num_fig)
file_out = [file 'row' '.png'];
plot(1:N,A(N/2,:));
axis([1 N 0 1.2]);
title([caption ' - The ' int2str(N/2) 'th row']);
print('-djpeg',file_out);

