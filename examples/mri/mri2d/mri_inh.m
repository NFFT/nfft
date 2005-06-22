N=256;   % points per row / column
M=75520; % number of knots

%load output_phantom_nfft.dat

%M=size(output_phantom_nfft,1)


%precompute_weights('output_phantom_nfft.dat',M);

disp(2);
for iter=1:1,
system(['./reconstruct_data_special_invers ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 1 1'])

load output_real.dat
figure(1);
imagesc(reshape(output_real,N,N));  
colormap(flipud(gray(256)));
colorbar;
%print('-dpng',['bild_iter=' int2str(iter) '.png'])

end
