f=phantom(128);
save 'input_data.dat' f -ascii

figure(1);
imagesc(f);
axis image
title('phantom');

!radon 128 256 256 5

Rf=load('sinogram_data.dat');

figure(2);
imagesc(Rf);
axis image
title('sinogram');
xlabel('offset');
ylabel('angle');


iRf=load('output_data.dat');

figure(3);
imagesc(iRf);
axis image
title('reconstructed image');

disp(sprintf('max(abs(f(:)-iRf(:))) = %g',max(abs(f(:)-iRf(:)))))