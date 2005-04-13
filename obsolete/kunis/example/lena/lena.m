% reconstruction of the lena image from one with lost data
% see http://www.nofiles.de/roots/lena/region.html

!make

N=256;
N_rec=256;
M=30000;

% take care of non periodic data set
Nc=5;

figure(1);
A=imread('lena256.jpg');
A(:,[1:Nc,end-Nc+1:end])=0;
A([1:Nc,end-Nc+1:end],:)=0;
imagesc(A);
colormap(gray(256));
axis([1+Nc,N-Nc,1+Nc,N-Nc]);
axis off

input_data=zeros(M,3);

[XI,YI]=meshgrid(1+Nc:N-Nc,1+Nc:N-Nc);
X=rand_subsample([XI(:),YI(:)],M);

figure(2);
B=zeros(N);
for j=1:M
  B(X(j,1),X(j,2))=A(X(j,1),X(j,2));
  input_data(j,3)=A(X(j,1),X(j,2));
end;
imagesc(B);
colormap(gray(256));
axis([1+Nc,N-Nc,1+Nc,N-Nc]);
axis off

input_data(:,1:2)=X/N-0.5;
save input_data.dat -ascii -double -tabs input_data

% use modulus dft as damping factors
input_w_hat=abs(ifftshift(fft2(fftshift(A))));
save input_w_hat.dat -ascii -double -tabs input_w_hat

% reconstruct phase
% change M,N here as well %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!lena 256 30000 30 > output_data.dat                           


load output_data.dat
f_hat=output_data(:,1)+i*output_data(:,2);

figure(3);
C=ifftshift(real(fft2(fftshift(reshape(f_hat,N_rec,N_rec)))));
C(:,[1:Nc,end-Nc+1:end])=0;
C([1:Nc,end-Nc+1:end],:)=0;
imagesc(uint8(C)');
colormap(gray(256));
axis([1+Nc,N-Nc,1+Nc,N-Nc]);
axis off
