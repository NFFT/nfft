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
%
% Visualizes the three dimensional phantom. Makes a pic of
% every plane and one plot of the N/2 row of the 10th plane.
function [] = visualize_data ( file, N, Z, num_fig, caption )

% load the real part of f
load output_real.dat

for z=0:Z-1,
  for y=0:N-1,
    for x =0:N-1,
      A(x+1,y+1,z+1)=output_real(z+1,y*N+ x+1);
     end
   end
end

% plot the three dimensional phantom
for j=0:Z-1,
  figure(2*num_fig-1)
  imagesc(A(:,:,j+1),[0 1]);
  colorbar;
  colormap(flipud(gray(256)));
  title(caption);
  if j<9,
    file_out =[file '0' int2str(j+1) '.jpg'];
  else
    file_out =[file int2str(j+1) '.jpg'];
  end
  print('-djpeg',file_out);
end

% plot the N/2 row of the 10 plane
figure(2*num_fig)
file_out = [file 'row' '.png'];
plot(1:N,A(:,N/2,Z/2));
axis([1 N 0 1.2]);
title([caption ' - The ' int2str(N/2) 'th row']);
print('-djpeg',file_out);
