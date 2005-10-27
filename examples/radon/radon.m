N=128;
T=2*N;
R=2*N;
it=5;

f=phantom(N);

fp = fopen('input_data.bin','wb+');
fwrite(fp,f','double');
fclose(fp);

figure(1);
imagesc(f);
axis image
title('phantom');

system(sprintf('radon %d %d %d %d',N,T,R,it));

fp = fopen('sinogram_data.bin','rb+');
Rf = fread(fp,[R,T],'double');
fclose(fp);

figure(2);
imagesc(Rf);
axis image
title('sinogram');
xlabel('angle');
ylabel('offset');


fp = fopen('output_data.bin','rb+');
iRf = fread(fp,[N,N],'double')';
fclose(fp);

figure(3);
imagesc(iRf);
axis image
title('reconstructed image');

disp(sprintf('max(abs(f(:)-iRf(:))) = %g',max(abs(f(:)-iRf(:)))))