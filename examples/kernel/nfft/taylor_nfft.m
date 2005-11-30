%% File: taylor_nfft.m
%%
%% Testing the nfft againt a Taylor expansion based version.
%%
%% Author Stefan Kunis
%%
%% References: [AnDa96], i.e., Chris Anderson and Marie Dahleh:
%% Rapid computation on the discrete Fourier transform
%%
%% Calls repeatedly the executable taylor_nfft.

to_pdf=0;

%%
%% Testing time vs. problem size.
%%
trials=10;
first=4;
last=20;
%system(sprintf('./taylor_nfft %d %d %d %d %f %f > taylor_nfft1.dat',0,first,last,trials,2,4));
data=load('taylor_nfft1.dat');

N=data(:,1);
avg_N=N(1:trials:end);
t_nfft=data(:,6);
avg_t_nfft=(mean(reshape(t_nfft,trials,last-first+1)))';
t_taylor=data(:,10);
avg_t_taylor=(mean(reshape(t_taylor,trials,last-first+1)))';

h=loglog(avg_N,avg_t_taylor,'ko',...
         avg_N,avg_t_nfft,'kx',...
         avg_N,10^-7*avg_N.*log(avg_N),'k:');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-6,10^-4,10^-2,1,10^2]);
set(gca,'XTick',[10^2,10^3,10^4,10^5,10^6]);
set(gca,'FontSize',20);
axis([N(1),N(end),10^-6,10^2]);

print temp.eps -deps
if(to_pdf)
  !ps2pdf temp.eps taylor_nfft1.pdf 
  !rm temp.eps
else
  !mv temp.eps taylor_nfft1.eps
end;

%%
%% Testing accuracy vs. cut-off/Taylor degree m.
%%
trials=20;
first=1;
last=20;
sigma_nfft=2;
sigma_taylor=2;
sigma_nfft1=1.5;
sigma_taylor1=1.5;
sigma_nfft2=16;
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

h=semilogy(first:last,max_e_taylor,'ko',...
           first:last,max_e_taylor_sigma,'k',...
           first:last,max_e_taylor_sigma2,'k--',...
           first:last,max_e_nfft,'kx',...
           first:last,max_e_nfft_sigma,'k-.',...
           first:last,max_e_nfft_sigma2,'k:');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-15,10^-10,10^-5,1]);
set(gca,'FontSize',20);
axis([first,last,10^-16,1]);

print temp.eps -deps
if(to_pdf)
  !ps2pdf temp.eps taylor_nfft2.pdf 
  !rm temp.eps
else
  !mv temp.eps taylor_nfft2.eps
end;

%%
%% Testing accuracy vs. time.
%%
t_nfft=data(:,6);
avg_t_nfft=(mean(reshape(t_nfft,trials,last-first+1)))';
t_taylor=data(:,10);
avg_t_taylor=(mean(reshape(t_taylor,trials,last-first+1)))';

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
