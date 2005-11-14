f=double(imread('phantom.png'));

figure(1);
imagesc(f);
axis image
title('phantom');

save 'input_data.dat' f -ascii

!polar_fft_test 128 256 256 5

f2=load('polar_fft_data.dat');

figure(2);
imagesc(real(f2));
axis image
title('reconstructed phantom (polar grid)');

!mpolar_fft_test 128 256 256 5

f2=load('mpolar_fft_data.dat');

figure(3);
imagesc(real(f2));
axis image
title('reconstructed phantom (modified polar grid)');

!linogram_fft_test 128 256 256 5

f2=load('linogram_fft_data.dat');

figure(4);
imagesc(real(f2));
axis image
title('reconstructed phantom (linogram grid)');
