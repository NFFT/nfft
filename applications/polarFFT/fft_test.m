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
N=64;
T=3*N; %T=5*N/2;
R=3*N/2;


f=fft2(phantom(N));
%f=imresize(im2double(imread('trui.png')),N/256);
fr=real(f); fi=imag(f);
save 'input_data_r.dat' fr -ascii -double
save 'input_data_i.dat' fi -ascii -double


if ispc
    cmd='polar_fft_test.exe';
else 
    cmd='./polar_fft_test';
end
system(sprintf('%s %d %d %d',cmd,N,T,R));

polar_fft_error = load('polar_fft_error.dat');
polar_ifft_error3 = load('polar_ifft_error3.dat');
polar_ifft_error6 = load('polar_ifft_error6.dat');
polar_ifft_error9 = load('polar_ifft_error9.dat');

figure(1);
h=semilogy(1:length(polar_fft_error),polar_fft_error);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22);
grid;
axis([0,10,10^-15,1]);
title('Test of the polar FFT');
xlabel('m');
ylabel('E_{max}');
print fig_polar_fft -deps2

figure(2);
it=0:10:100;
h=plot(it,polar_ifft_error3,'-',it,polar_ifft_error6,'--',it,polar_ifft_error9,'-.');
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22);
grid;
title('Test of the inverse polar FFT');
xlabel('iterations');
ylabel('E_{max}');
legend('m=3','m=6','m=9');
print fig_ipolar_fft -deps2

disp(sprintf('\n'));


if ispc
    cmd_mpolar_fft_test='mpolar_fft_test.exe';
else 
    cmd_mpolar_fft_test='./mpolar_fft_test';
end
system(sprintf('%s %d %d %d',cmd_mpolar_fft_test,N,T,R));

mpolar_fft_error = load('mpolar_fft_error.dat');
mpolar_ifft_error3 = load('mpolar_ifft_error3.dat');
mpolar_ifft_error6 = load('mpolar_ifft_error6.dat');
mpolar_ifft_error9 = load('mpolar_ifft_error9.dat');

figure(3);
h=semilogy(1:length(mpolar_fft_error),mpolar_fft_error);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22);
grid;
axis([0,10,10^-15,1]);
title('Test of the mpolar FFT');
xlabel('m');
ylabel('E_{max}');
print fig_mpolar_fft -deps2

figure(4);
it=0:2:20;
h=semilogy(it,mpolar_ifft_error3,'-',it,mpolar_ifft_error6,'--',it,mpolar_ifft_error9,'-.');
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22);
grid;
title('Test of the inverse mpolar FFT');
xlabel('iterations');
ylabel('E_{max}');
legend('m=3','m=6','m=9');
print fig_impolar_fft -deps2

disp(sprintf('\n'));


if ispc
    cmd_linogram_fft_test='linogram_fft_test.exe';
else 
    cmd_linogram_fft_test='./linogram_fft_test';
end
system(sprintf('%s %d %d %d',cmd_linogram_fft_test,N,T,R));

linogram_fft_error = load('linogram_fft_error.dat');
linogram_ifft_error3 = load('linogram_ifft_error3.dat');
linogram_ifft_error6 = load('linogram_ifft_error6.dat');
linogram_ifft_error9 = load('linogram_ifft_error9.dat');

figure(5);
h=semilogy(1:length(linogram_fft_error),linogram_fft_error);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22);
grid;
axis([0,10,10^-15,1]);
title('Test of the linogram FFT');
xlabel('m');
ylabel('E_{max}');
print fig_lino_fft -deps2

figure(6);
it=0:2:20;
h=semilogy(it,linogram_ifft_error3,'-',it,linogram_ifft_error6,'--',it,linogram_ifft_error9,'-.');
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22);
grid;
title('Test of the inverse linogram FFT');
xlabel('iterations');
ylabel('E_{max}');
legend('m=3','m=6','m=9');
print fig_ilino_fft -deps2
