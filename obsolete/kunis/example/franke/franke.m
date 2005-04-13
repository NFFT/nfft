!make

N=256;
M=10000;

FrankeRange=2;
NfftRange=1/4;
ratio=FrankeRange/NfftRange;

X=2*NfftRange*rand(M,1)-NfftRange;
Y=2*NfftRange*rand(M,1)-NfftRange;
Z=f1(ratio*X,ratio*Y);

input_data=[X,Y,Z];

save input_data.dat -ascii -double -tabs input_data

!franke 256 10000 20 > output_data.dat

x_rec=-1/2:1/N:1/2-1/N;
[X_rec,Y_rec]=meshgrid(x_rec,x_rec);
Z_rec=f1(ratio*X_rec,ratio*Y_rec);

load output_data.dat
output=output_data(:,1)+i*output_data(:,2);
A=reshape(output,N,N);

Z_rec_nfft=real(ifftshift(fft2(fftshift(A))));

figure(1); surf(X_rec,Y_rec,Z_rec); shading interp; colorbar; title('original function'); print franke1.eps -depsc
figure(2); contourf(X_rec,Y_rec,Z_rec_nfft); axis([-1/4,1/4,-1/4,1/4]); colorbar; title('reconstructed function'); print franke2.eps -depsc
figure(3); contourf(X_rec,Y_rec,log(abs(Z_rec-Z_rec_nfft))); axis([-1/4,1/4,-1/4,1/4]); colorbar; title('log of absolute residual'); print franke3.eps -depsc