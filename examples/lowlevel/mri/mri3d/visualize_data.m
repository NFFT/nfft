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
title([caption ' - The 64th row']);
print('-djpeg',file_out);
