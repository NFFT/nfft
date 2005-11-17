N=128;
f=phantom(N);
%grid='polar'; T=2.5*N; R=1.5*N; it=5;
grid='linogram'; T=2*N; R=2*N; it=5;

fp = fopen('input_data.bin','wb+');
fwrite(fp,f','double');
fclose(fp);

figure(1);
imagesc(f);
axis image
title('phantom');

system(sprintf('./radon %s %d %d %d %d',grid,N,T,R,it));

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

disp(sprintf('max(abs(f(:)-iRf(:))) = %e',max(abs(f(:)-iRf(:)))))
