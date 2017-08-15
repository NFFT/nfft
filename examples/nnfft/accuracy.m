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
m=0:14;
load kaiser_bessel.dat
load sinc_power.dat
load b_spline.dat
load gaussian.dat

figure(1);
h=semilogy(m,kaiser_bessel(m+1,1),'k',m,kaiser_bessel(m+1,1),'ko',...
           m,sinc_power(m+1,1),'k',m,sinc_power(m+1,1),'kx',...
           m,b_spline(m+1,1),'k',m,b_spline(m+1,1),'k+',...
           m,gaussian(m+1,1),'k',m,gaussian(m+1,1),'k^');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print accuracy1.eps -deps

figure(2);
h=semilogy(m,kaiser_bessel(m+1,2),'k',m,kaiser_bessel(m+1,2),'ko',...
           m,sinc_power(m+1,2),'k',m,sinc_power(m+1,2),'kx',...
           m,b_spline(m+1,2),'k',m,b_spline(m+1,2),'k+',...
           m,gaussian(m+1,2),'k',m,gaussian(m+1,2),'k^');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print accuracy2.eps -deps

figure(3);
h=semilogy(m,kaiser_bessel(m+16,1),'k',m,kaiser_bessel(m+16,1),'ko',...
           m,sinc_power(m+16,1),'k',m,sinc_power(m+16,1),'kx',...
           m,b_spline(m+16,1),'k',m,b_spline(m+16,1),'k+',...
           m,gaussian(m+16,1),'k',m,gaussian(m+16,1),'k^');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print accuracy3.eps -deps

figure(4);
h=semilogy(m,kaiser_bessel(m+16,2),'k',m,kaiser_bessel(m+16,2),'ko',...
           m,sinc_power(m+16,2),'k',m,sinc_power(m+16,2),'kx',...
           m,b_spline(m+16,2),'k',m,b_spline(m+16,2),'k+',...
           m,gaussian(m+16,2),'k',m,gaussian(m+16,2),'k^');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print accuracy4.eps -deps

figure(5);
h=semilogy(m,kaiser_bessel(m+31,1),'k',m,kaiser_bessel(m+31,1),'ko',...
           m,sinc_power(m+31,1),'k',m,sinc_power(m+31,1),'kx',...
           m,b_spline(m+31,1),'k',m,b_spline(m+31,1),'k+',...
           m,gaussian(m+31,1),'k',m,gaussian(m+31,1),'k^');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print accuracy5.eps -deps

figure(6);
h=semilogy(m,kaiser_bessel(m+31,2),'k',m,kaiser_bessel(m+31,2),'ko',...
           m,sinc_power(m+31,2),'k',m,sinc_power(m+31,2),'kx',...
           m,b_spline(m+31,2),'k',m,b_spline(m+31,2),'k+',...
           m,gaussian(m+31,2),'k',m,gaussian(m+31,2),'k^');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print accuracy6.eps -deps
