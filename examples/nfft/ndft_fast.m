%% File: ndft_fast.m
%%
%% Testing ndft, Horner-like ndft, fully precomputed ndft and nfft.
%% Testing time vs. problem size.
%%
%% Author Stefan Kunis
%%
%% References: Time and memory requirements of the Nonequispaced FFT
%%
%% Calls repeatedly the executable ndft_fast.

to_pdf=0;

trials=10;
first=4;
last=17;
%system(sprintf('./ndft_fast %d %d %d %d > ndft_fast.data',0,first,last,trials));
data=load('ndft_fast.data');

N=data(1:trials:end,1);
t_ndft=(mean(reshape(data(:,3),trials,last-first+1)))';
t_horner=(mean(reshape(data(:,4),trials,last-first+1)))';
t_pre_full=(mean(reshape(data(:,5),trials,last-first+1)))';
t_nfft=(mean(reshape(data(:,6),trials,last-first+1)))';

h=loglog(N,t_ndft,'k',...
         N,t_horner,'k--',...
         N,t_pre_full,'k-.',...
         N,10^-8*N.^2,'k:',...
         N,t_nfft,'kd');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-6,10^-4,10^-2,1,10^2]);
set(gca,'XTick',[10^2,10^3,10^4,10^5,10^6]);
set(gca,'FontSize',20);
axis([N(1),N(end),10^-6,10^2]);

print temp.eps -deps
if(to_pdf)
  !ps2pdf temp.eps ndft_fast.pdf 
  !rm temp.eps
else
  !mv temp.eps ndft_fast.eps
end;
