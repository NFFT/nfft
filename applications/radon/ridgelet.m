% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%
%file: ridgelet.m
%
%  Simple test program for denoising an image
%  by hard thresholding the ridgelet coefficients.
%  It uses the NFFT-based discrete Radon transform
%  and translationinvariant discrete Wavelet transform.
%  Requires Matlab-Toolbox WaveLab802.
%
%references: D. Donoho et al.,
%  WaveLab802 - A collection of Matlab functions.
%  http://www-stat.stanford.edu/~wavelab/, 1999.
%
%author: Markus Fenn
%date: April 2006

N=128;
grid='linogram'; T=2*N; R=2*N; it=5;  %grid and parameters for Radon transform
SNR=1;                                %noise level
qmf=MakeONFilter('Battle',3);         %filter for translation invariant DWT
threshold=17;                         %threshhold for denoising

%original input image
f=ones(N,N);
for j=1:N
  for k=j:(N-j)
    f(j,k)=0.0;
  end
end
f((5*N/8):(7*N/8),(N/4):(3*N/4))=0.0;

%plot original image
figure(1);
colormap(gray);
imagesc(f);
axis image
title('original image');

%add gaussian white noise of given SNR
noise=randn(size(f));
noise_factor = norm(f - mean(f(:)),'fro')/10^(SNR/20);
f2 = f + noise_factor * noise/norm(noise,'fro');

%write noisy image to file
fp = fopen('input_data.bin','wb+');
fwrite(fp,f2','double');
fclose(fp);

%plot noisy image
figure(2);
colormap(gray);
imagesc(f2);
axis image
title(sprintf('noisy image (SNR=%g)',SNR));

%compute Radon transform by C-program
if ispc
    cmd='radon.exe';
else 
    cmd='./radon';
end
system(sprintf('%s %s %d %d %d %d',cmd,grid,N,T,R));

%read result from file
fp = fopen('sinogram_data.bin','rb+');
Rf2 = fread(fp,[R,T],'double');
fclose(fp);

%compute translation invariant DWT for every direction
Rf2d=zeros(size(Rf2));
for k=1:T
  TItable = FWT_TI(Rf2(:,k)', 4, qmf);
  [nrow,ncol] = size(TItable);
  for j = 2:ncol
      TItable(:,j) = TItable(:,j) .* (abs(TItable(:,j)) > threshold);
  end
  temp = IWT_TI(TItable,qmf);
  Rf2d(:,k)=temp';
end

%write result to file
fp = fopen('sinogram_data.bin','wb+');
fwrite(fp,Rf2d,'double');
fclose(fp);

%compute inverse Radon transform by C-program
if ispc
    cmd='inverse_radon.exe';
else 
    cmd='./inverse_radon';
end
system(sprintf('%s %s %d %d %d %d',cmd,grid,N,T,R,it));

%read result
fp = fopen('output_data.bin','rb+');
iRf2d = fread(fp,[N,N],'double')';
fclose(fp);

figure(3);
imagesc(iRf2d);
colormap(gray);
axis image
title(sprintf('reconstructed image (SNR=%g)',20*log10( norm(f - mean(f(:)),'fro') / norm(f-iRf2d,'fro') )));

%end: ridgelet.m
