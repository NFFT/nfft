%% File: flags.m
%%
%% Testing the Testing different precomputation schemes for the nfft.
%%
%% Author Stefan Kunis
%%
%% References: [flags], i.e., Draft:
%% Time and memory requirements of the Nonequispaced FFT
%%
%% Calls repeatedly the executable flags.

%%
%% Testing time vs. problem size.
%%
trials=2;
first=4;
last=20;
system(sprintf('./flags %d %d %d %d %f %f > flags1.dat',0,first,last,trials,1,4));
data=load('flags1.dat');

N=data(:,2);
avg_N=N(1:trials:end);
t_D=data(:,5);
avg_t_D=(mean(reshape(t_D,trials,last-first+1)))';
t_pre_phi_hut=data(:,6);
avg_t_pre_phi_hut=(mean(reshape(t_pre_phi_hut,trials,last-first+1)))';
t_fftw=data(:,7);
avg_t_fftw=(mean(reshape(t_fftw,trials,last-first+1)))';
t_B=data(:,8);
avg_t_B=(mean(reshape(t_B,trials,last-first+1)))';
t_fg_psi=data(:,9);
avg_t_fg_psi=(mean(reshape(t_fg_psi,trials,last-first+1)))';
t_pre_lin_psi=data(:,10);
avg_t_pre_lin_psi=(mean(reshape(t_pre_lin_psi,trials,last-first+1)))';
t_pre_fg_psi=data(:,11);
avg_t_pre_fg_psi=(mean(reshape(t_pre_fg_psi,trials,last-first+1)))';
t_pre_psi=data(:,12);
avg_t_pre_psi=(mean(reshape(t_pre_psi,trials,last-first+1)))';
t_pre_full_psi=data(:,5);
avg_t_pre_full_psi=(mean(reshape(t_pre_full_psi,trials,last-first+1)))';

h=loglog(N,t_D,'go',avg_N,avg_t_D,'g',...
         N,t_pre_phi_hut,'go',avg_N,avg_t_pre_phi_hut,'g:',...
         N,t_fftw,'ro',avg_N,avg_t_fftw,'r',...
         N,t_B,'ko',avg_N,avg_t_B,'k',...
         N,t_fg_psi,'kx',avg_N,avg_t_fg_psi,'k:',...
         N,t_pre_lin_psi,'k+',avg_N,avg_t_pre_lin_psi,'k-.',...
         N,t_pre_fg_psi,'k*',avg_N,avg_t_pre_fg_psi,'k--',...
         N,t_pre_psi,'ks',avg_N,avg_t_pre_psi,'k:',...
         N,t_pre_full_psi,'k*',avg_N,avg_t_pre_full_psi,'k-.');
  
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-6,10^-4,10^-2,1,10^2]);
set(gca,'XTick',[10^2,10^3,10^4,10^5,10^6]);
set(gca,'FontSize',20);
axis([N(1),N(end),10^-6,10^0]);

print temp.eps -deps
!ps2pdf temp.eps flags1.pdf 
!rm temp.eps

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
system(sprintf('./taylor_nfft %d %d %d %d %f %f > taylor_nfft2a.dat',1,first,...
	       last,trials,sigma_nfft,sigma_taylor));
data=load('taylor_nfft2a.dat');

m_nfft=data(:,5);
e_nfft=data(:,7);
max_e_nfft=(max(reshape(e_nfft,trials,last-first+1)))';

m_taylor=data(:,9);
e_taylor=data(:,11);
max_e_taylor=(max(reshape(e_taylor,trials,last-first+1)))';

% small sigma
system(sprintf('./taylor_nfft %d %d %d %d %f %f > taylor_nfft2b.dat',1,first,...
	       last,trials,sigma_nfft1,sigma_taylor1));
data1=load('taylor_nfft2b.dat');
max_e_nfft_sigma=(max(reshape(data1(:,7),trials,last-first+1)))';
max_e_taylor_sigma=(max(reshape(data1(:,11),trials,last-first+1)))';

% large sigma
system(sprintf('./taylor_nfft %d %d %d %d %f %f > taylor_nfft2c.dat',1,first,...
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

