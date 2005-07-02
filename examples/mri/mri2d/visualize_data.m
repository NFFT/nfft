function [] = visualize_data ( file, N )

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
imagesc(A,[0 1.2]);
colormap((gray(256)));
colorbar;
file_out =[file '.png'];
print('-dpng',file_out);

% plot the N/2 row
file_out = [file 'row' '.png'];
plot(1:N,A(N/2,:));
axis([1 N 0 1.2]);
print('-dpng',file_out);

