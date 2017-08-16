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
%file: radon.m
%
%  Simple test program for NFFT-based discrete
%  Radon transform and its inverse.
%
%  See radon.c, inverse_radon.c and refman.pdf for details.
%
%author: Markus Fenn
%date: February 2006

N=128;
f=phantom(N);
%grid='polar'; T=2.5*N; R=1.5*N; it=5;
grid='linogram'; T=2*N; R=2*N; it=5;

%plot input image
figure(1);
imagesc(f);
axis image
title('phantom');

%write input to file
fp = fopen('input_data.bin','wb+');
fwrite(fp,f','double');
fclose(fp);

%compute Radon transform by C-program
if ispc
    cmd='radon.exe';
else 
    cmd='./radon';
end
system(sprintf('%s %s %d %d %d',cmd,grid,N,T,R));

%read result from file
fp = fopen('sinogram_data.bin','rb+');
Rf = fread(fp,[R,T],'double');
fclose(fp);

%plot sinogram
figure(2);
imagesc(Rf);
axis image
title('sinogram');
xlabel('angle');
ylabel('offset');

%compute inverse Radon transform by C-program
if ispc
    cmd='inverse_radon.exe';
else 
    cmd='./inverse_radon';
end
system(sprintf('%s %s %d %d %d %d',cmd,grid,N,T,R,it));

%read result from file
fp = fopen('output_data.bin','rb+');
iRf = fread(fp,[N,N],'double')';
fclose(fp);

%plot reconstructed image
figure(3);
imagesc(iRf);
axis image
title('reconstructed image');

%compute error
disp(sprintf('max(abs(f(:)-iRf(:))) = %e',max(abs(f(:)-iRf(:)))))

%end: radon.m
