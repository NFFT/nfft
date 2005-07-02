N=256;   % points per row / column
arms=32;
M=construct_knots_spiral(256,arms);

%construct_readout_time( M, 2, arms, 0.00402542373 );
construct_readout_time( M, 2, arms, 0.004 );

construct_inh(N);

out=reshape(phantom(N),1,N*N);
save input_f.dat -ascii out 

system(['./construct_data_inhomogen ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)])


precompute_weights('output_phantom_nfft.dat',M);

system(['./reconstruct_data_special_invers ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 10 1'])


load output_real.dat
load output_imag.dat
output_real=abs(output_real+i*output_imag);

figure(1);
imagesc(reshape(output_real,N,N));  
colormap(flipud(gray(256)));
colorbar;
%print('-dpng',['bild_iter=' int2str(iter) '.png'])

system(['./reconstruct_data_2d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 4 1'])

load output_real.dat
load output_imag.dat
output_real=abs(output_real+i*output_imag);

figure(2);
imagesc(reshape(output_real,N,N));
colormap(flipud(gray(256)));
colorbar;
%print('-dpng',['bild_iter=' int2str(iter) '.png'])

