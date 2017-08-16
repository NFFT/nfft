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
output_error

a=real(delta);
C=exp(-a*pi^2/abs(delta)^2);

figure(1);
h=semilogy( error(:,1), error(:,2),'k-',...
  	  error(:,1), error(:,3),'k--',...
	  error(:,1), error(:,4),'k-.',...
error(:,1), 2*exp(-a/4)*(1+1/a) + sqrt(pi/abs(delta)) * exp(-pi^2*a*error(:,1).^2/(4*abs(delta)^2)).*(1+abs(delta)^2./(error(:,1)*pi^2*a)),'k:');
set(h,'LineWidth',1.8); set(h,'Markersize',10); set(gca,'FontSize',25);
axis([min(error(:,1)),max(error(:,1)),10^-18,1]);
print gauss1.eps -deps

figure(2);
h=semilogy( error(:,1), error(:,5),'k-',...
	  error(:,1), error(:,6),'k--',...
	  error(:,1), error(:,7),'k-.',...
	  error(:,1), 8*abs(delta)^2./(pi^2*a*error(:,1)).*exp(-a/4).*(1+1./error(:,1)) + sqrt(pi/abs(delta)) * exp(-pi^2*a*error(:,1).^2/(4*abs(delta)^2)).*(1+2*abs(delta)^2./(error(:,1)*pi^2*a)),'k:');
set(h,'LineWidth',1.8); set(h,'Markersize',10); set(gca,'FontSize',25);
axis([min(error(:,1)),max(error(:,1)),10^-18,1]);
print gauss2.eps -deps



output_error_p

figure(3);
h=semilogy( error(:,1), error(:,2),'k-',...
	    error(:,1), error(:,3),'k--',...
	    error(:,1), error(:,4),'k-.');
set(h,'LineWidth',1.8); set(h,'Markersize',10); set(gca,'FontSize',25);
axis([min(error(:,1)),max(error(:,1)),10^-18,1]);
print gauss6.eps -deps
