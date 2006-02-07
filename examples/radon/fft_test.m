N=128;
T=2*N;
R=N;
it=10;

f=phantom(N);

figure(1);
subplot(2,2,1);
imagesc(f);
axis image
title('phantom');

save 'input_data.dat' f -ascii

system(sprintf('./polar_fft_test %d %d %d %d',N,T,R,it));

f2=load('polar_fft_data.dat');

subplot(2,2,2);
imagesc(real(f2));
axis image
title('reconstructed phantom (polar grid)');

system(sprintf('./mpolar_fft_test %d %d %d %d',N,T,R,it));

f2=load('mpolar_fft_data.dat');

subplot(2,2,3);
imagesc(real(f2));
axis image
title('reconstructed phantom (modified polar grid)');

system(sprintf('./linogram_fft_test %d %d %d %d',N,T,R,it));

f2=load('linogram_fft_data.dat');

subplot(2,2,4);
imagesc(real(f2));
axis image
title('reconstructed phantom (linogram grid)');
