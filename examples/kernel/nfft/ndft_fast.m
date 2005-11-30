%% File: ndft_fast.m
%%
%% Testing ndft, Horner-like ndft, fully precomputed ndft and nfft.
%%
%% Author Stefan Kunis
%%
%% References: [BaMi], i.e., Bagchi Mitra XXXXXXXX:
%%
%% Calls repeatedly the executable ndft_fast.

to_pdf=0;

%%
%% Testing time vs. problem size.
%%
trials=2;%10;
first=4;
last=16;
%system(sprintf('./ndft_fast %d %d %d %d > ndft_fast1.dat',0,first,last,trials));
data=load('ndft_fast1.dat');

N=data(:,1);
avg_N=N(1:trials:end);
t_ndft=data(:,3);
avg_t_ndft=(mean(reshape(t_ndft,trials,last-first+1)))';
t_horner=data(:,4);
avg_t_horner=(mean(reshape(t_horner,trials,last-first+1)))';
t_pre_full=data(:,5);
avg_t_pre_full=(mean(reshape(t_pre_full,trials,last-first+1)))';
t_nfft=data(:,6);
avg_t_nfft=(mean(reshape(t_nfft,trials,last-first+1)))';

h=loglog(avg_N,avg_t_ndft,'ko',...
         avg_N,avg_t_horner,'kx',...
         avg_N,avg_t_pre_full,'k+',...
         avg_N,10^-8*avg_N.^2,'k:',...
         avg_N,avg_t_nfft,'kd');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-6,10^-4,10^-2,1,10^2]);
set(gca,'XTick',[10^2,10^3,10^4,10^5,10^6]);
set(gca,'FontSize',20);
axis([N(1),N(end),10^-6,10^2]);

print temp.eps -deps
if(to_pdf)
  !ps2pdf temp.eps ndft_fast1.pdf 
  !rm temp.eps
else
  !mv temp.eps ndft_fast1.eps
end;
