function [] = visualize_phantom ( file, N )

load output_real.dat

A=reshape(output_real,N,N);

% plot the two dimensional phantom
figure(1)
imagesc(A,[0 1]);
colormap(flipud(gray(256)));
colorbar;
file_out =[file '.png'];
print('-dpng',file_out);

% plot the N/2 row
file_out = [file 'row' '.png'];
plot(1:N,A(N/2,:));
axis([1 N 0 1]);
print('-dpng',file_out);

