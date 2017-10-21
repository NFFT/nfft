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
%% File: flags1.m
%%
%% Testing different precomputation schemes for the nfft.
%% 1. Computation time vs. problem size
%%
%% Author Stefan Kunis
%%
%% References: Time and memory requirements of the Nonequispaced FFT
%%
%% Loads flags.data0.gaussian, ...

% see flags.readme
trials=10;
first=4;
last=20;

% load already computed data, times for precomputed D/B differ only negligible
dg=load('flags.data0.gaussian');
dk=load('flags.data0.kaiser_bessel');
ds=load('flags.data0.sinc_power');
db=load('flags.data0.b_spline');

N=dg(1:trials:end,2);

t_Dg=(mean(reshape(dg(:,6),trials,last-first+1)))';
t_Dk=(mean(reshape(dk(:,6),trials,last-first+1)))';
t_Ds=(mean(reshape(ds(:,6),trials,last-first+1)))';
t_Db=(mean(reshape(db(:,6),trials,last-first+1)))';
t_pre_phi_hut=(mean(reshape(dg(:,7),trials,last-first+1)))';

t_fftw=(mean(reshape(dg(:,8),trials,last-first+1)))';

t_Bg=(mean(reshape(dg(:,9),trials,last-first+1)))';
t_Bk=(mean(reshape(dk(:,9),trials,last-first+1)))';
t_Bs=(mean(reshape(ds(:,9),trials,last-first+1)))';
t_Bb=(mean(reshape(db(:,9),trials,last-first+1)))';

t_fg_psi=(mean(reshape(dg(:,10),trials,last-first+1)))';
t_pre_lin_psi=(mean(reshape(dg(:,11),trials,last-first+1)))';
t_pre_fg_psi=(mean(reshape(dg(:,12),trials,last-first+1)))';
t_pre_psi=(mean(reshape(dg(:,13),trials,last-first+1)))';
t_pre_full_psi=(mean(reshape(dg(:,14),trials,last-first+1)))';

h=loglog(N,t_Dg,'k',...
         N,t_Dk,'k--',...
         N,t_Ds,'k-.',...
         N,t_Db,'k:',...
         N,t_pre_phi_hut,'k.');
  
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-6,10^-4,10^-2,1,10^2]);
set(gca,'XTick',[10^2,10^3,10^4,10^5,10^6]);
set(gca,'FontSize',20);
axis([N(1),N(end),10^-6,10^1]);

if(to_pdf)
  !ps2pdf temp.eps flags1.pdf 
  !rm temp.eps
else
  !mv temp.eps flags1.eps
end;

return

%%
%% Testing accuracy vs. cut-off/Taylor degree m.
%%
trials=20;
first=1;
last=20;
sigma_nfft=2;
sigma_taylor=4;
sigma_nfft1=1.25;
sigma_taylor1=1.5;
sigma_nfft2=8;
sigma_taylor2=16;

% typical sigma
if ispc
    cmd='taylor_nfft.exe';
else 
    cmd='./taylor_nfft';
end
system(sprintf('%s %d %d %d %d %f %f > taylor_nfft2a.dat',cmd,1,first,...
	       last,trials,sigma_nfft,sigma_taylor));
data=load('taylor_nfft2a.dat');

m_nfft=data(:,5);
e_nfft=data(:,7);
max_e_nfft=(max(reshape(e_nfft,trials,last-first+1)))';

m_taylor=data(:,9);
e_taylor=data(:,11);
max_e_taylor=(max(reshape(e_taylor,trials,last-first+1)))';

% small sigma
system(sprintf('%s %d %d %d %d %f %f > taylor_nfft2b.dat',cmd,1,first,...
	       last,trials,sigma_nfft1,sigma_taylor1));
data1=load('taylor_nfft2b.dat');
max_e_nfft_sigma=(max(reshape(data1(:,7),trials,last-first+1)))';
max_e_taylor_sigma=(max(reshape(data1(:,11),trials,last-first+1)))';

% large sigma
system(sprintf('%s %d %d %d %d %f %f > taylor_nfft2c.dat',cmd,1,first,...
	       last,trials,sigma_nfft2,sigma_taylor2));
data2=load('taylor_nfft2c.dat');
max_e_nfft_sigma2=(max(reshape(data2(:,7),trials,last-first+1)))';
max_e_taylor_sigma2=(max(reshape(data2(:,11),trials,last-first+1)))';

h=semilogy(m_nfft,e_nfft,'k+',first:last,max_e_nfft,'k--',first:last,...
	   max_e_nfft_sigma,'k:',first:last,max_e_nfft_sigma2,'k:',...
	   m_taylor,e_taylor,'kx',first:last,max_e_taylor,'k',first:last,...
	   max_e_taylor_sigma,'k-.',first:last,max_e_taylor_sigma2,'k-.');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-15,10^-10,10^-5,1]);
set(gca,'FontSize',20);
axis([first,last,10^-16,1]);

print temp.eps -deps
!ps2pdf temp.eps taylor_nfft2.pdf 
!rm temp.eps


%%
%% Testing accuracy vs. time.
%%
t_nfft=data(:,6);
t_taylor=data(:,10);

h=semilogy(t_nfft,e_nfft,'k+',t_taylor,e_taylor,'kx');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-15,10^-10,10^-5,1]);
set(gca,'FontSize',20);
axis([min([t_nfft;t_taylor]),max([t_nfft;t_taylor]),10^-16,1]);

print temp.eps -deps
!ps2pdf temp.eps taylor_nfft3.pdf 
!rm temp.eps

