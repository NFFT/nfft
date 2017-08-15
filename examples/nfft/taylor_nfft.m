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

%% File: taylor_nfft.m
%%
%% Testing the nfft againt a Taylor expansion based version.
%%
%% Author Stefan Kunis
%%
%% References: Time and memory requirements of the Nonequispaced FFT
%%
%% Calls repeatedly the executable taylor_nfft.

to_pdf=0;

%%
%% Testing time vs. problem size.
%%
trials=10;
first=4;
last=22;
if ispc
    cmd='taylor_nfft.exe';
else 
    cmd='./taylor_nfft';
end
system(sprintf('%s %d %d %d %d %f %f > taylor_nfft.data0',cmd,0,first,last,trials,2,4));
data=load('taylor_nfft.data0');

N=data(1:trials:end,1);
t_nfft=(mean(reshape(data(:,6),trials,last-first+1)))';
t_taylor=(mean(reshape(data(:,10),trials,last-first+1)))';

h=loglog(N,t_taylor,'k',...
         N,t_nfft,'k--',...
         N,10^-7*N.*log(N),'k:');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-6,10^-4,10^-2,1,10^2]);
set(gca,'XTick',[10^2,10^3,10^4,10^5,10^6]);
set(gca,'FontSize',20);
axis([N(1),N(end),10^-6,10^2]);

print temp.eps -deps
if(to_pdf)
  !ps2pdf temp.eps taylor_nfft0.pdf 
  !rm temp.eps
else
  !mv temp.eps taylor_nfft0.eps
end;

return
%%
%% Testing accuracy vs. cut-off/Taylor degree m.
%%
trials=10;
first=1;
last=20;
sigma_nfft_a=2;
sigma_taylor_a=2;
sigma_nfft_b=1.5;
sigma_taylor_b=1.5;
sigma_nfft_c=16;
sigma_taylor_c=16;

% typical sigma
system(sprintf('%s %d %d %d %d %f %f > taylor_nfft.data1a',cmd,1,first,...
	       last,trials,sigma_nfft_a,sigma_taylor_a));
data=load('taylor_nfft.data1a');

m_nfft=data(1:trials:end,5);
e_nfft_a=(max(reshape(data(:,7),trials,last-first+1)))';
t_nfft_a=(mean(reshape(data(:,6),trials,last-first+1)))';

m_taylor=data(1:trials:end,9);
e_taylor_a=(max(reshape(data(:,11),trials,last-first+1)))';
t_taylor_a=(mean(reshape(data(:,10),trials,last-first+1)))';

% small sigma
system(sprintf('%s %d %d %d %d %f %f > taylor_nfft.data1b',cmd,1,first,...
	       last,trials,sigma_nfft_b,sigma_taylor_b));
data=load('taylor_nfft.data1b');

e_nfft_b=(max(reshape(data1(:,7),trials,last-first+1)))';
e_taylor_b=(max(reshape(data1(:,11),trials,last-first+1)))';

% large sigma
system(sprintf('%s %d %d %d %d %f %f > taylor_nfft.data1c',cmd,1,first,...
	       last,trials,sigma_nfft_c,sigma_taylor_c));
data=load('taylor_nfft.data1c');

e_nfft_c=(max(reshape(data2(:,7),trials,last-first+1)))';
e_taylor_c=(max(reshape(data2(:,11),trials,last-first+1)))';

h=semilogy(m_taylor,e_taylor_a,'kd',...
	   m_nfft,e_nfft_a,'kx',...
	   m_taylor,e_taylor_b,'k',...
	   m_nfft,e_nfft_b,'k--',...
	   m_taylor,e_taylor,'k-.',...
	   m_nfft,e_nfft,'k:');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-15,10^-10,10^-5,1]);
set(gca,'FontSize',20);
axis([m(1),m(end),10^-13,10^-1]);

print temp.eps -deps
if(to_pdf)
  !ps2pdf temp.eps taylor_nfft1.pdf 
  !rm temp.eps
else
  !mv temp.eps taylor_nfft1.eps
end;

%%
%% Testing accuracy vs. time.
%%


h=semilogy(avg_t_taylor,max_e_taylor,'ko',...
           avg_t_nfft,max_e_nfft,'kx');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-15,10^-10,10^-5,1]);
set(gca,'FontSize',20);
axis([min([t_nfft;t_taylor]),max([t_nfft;t_taylor]),10^-16,1]);

print temp.eps -deps
if(to_pdf)
  !ps2pdf temp.eps taylor_nfft3.pdf 
  !rm temp.eps
else
  !mv temp.eps taylor_nfft3.eps
end;
