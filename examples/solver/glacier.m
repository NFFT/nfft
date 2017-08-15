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
N=256;
border_eps=0.1;

M=8345;

load vol87.dat -ascii
input_data=vol87(2:end,:);

x_range=max(input_data(:,1))-min(input_data(:,1));
input_data(:,1)=(input_data(:,1)-min(input_data(:,1))) /x_range*(1-2*border_eps)-0.5+border_eps;

y_range=max(input_data(:,2))-min(input_data(:,2));
input_data(:,2)=(input_data(:,2)-min(input_data(:,2))) /y_range*(1-2*border_eps)-0.5+border_eps;

save input_data.dat -ascii -double -tabs input_data

if ispc
    cmd='glacier.exe';
else 
    cmd='./glacier';
end
system(sprintf('%s %d %d > output_data.dat',cmd,N,M));

load output_data.dat
f_hat=output_data(:,1)+i*output_data(:,2);

figure(1);
C=ifftshift(real(fft2(fftshift(reshape(f_hat,N,N)))));
x=-0.5:1/N:0.5-1/N;
[X,Y]=meshgrid(x,x);
Border=ceil(border_eps*N);
X=X(Border+1:N-Border,Border+1:N-Border);
Y=Y(Border+1:N-Border,Border+1:N-Border);
C=C(Border+1:N-Border,Border+1:N-Border);
surfl(X,Y,C);
axis([-0.5+border_eps,0.5-border_eps,-0.5+border_eps,0.5-border_eps,min(input_data(:,3)),max(C(:))]);
colormap(bone);
shading interp
axis off
view([30,30])
print glacier1.eps -deps

figure(2);
contour(X,Y,C,sort(input_data(:,3)));
axis([-0.5+border_eps,0.5-border_eps,-0.5+border_eps,0.5-border_eps]);
axis off
hold on
plot(input_data(:,1),input_data(:,2),'k.');
colormap(gray(256));
hold off
print glacier2.eps -deps
