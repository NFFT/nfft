%% File: ndft_fast.m
%%
%% Testing ndft, Horner-like ndft, and fully precomputed ndft.
%%
%% Author Stefan Kunis
%%
%% References: [BaMi], i.e., Bagchi Mitra XXXXXXXX:
%%
%% Calls repeatedly the executable ndft_fast.

%%
%% Testing time vs. problem size.
%%
trials=10;
first=1;
last=16;
system(sprintf('./ndft_fast %d %d %d %d > ndft_fast1.dat',0,first,last,trials));
data=load('ndft_fast1.dat');

N=data(:,1);
avg_N=N(1:trials:end);
t_ndft=data(:,3);
avg_t_ndft=(mean(reshape(t_ndft,trials,last-first+1)))';
t_horner=data(:,4);
avg_t_horner=(mean(reshape(t_horner,trials,last-first+1)))';
t_pre_full=data(:,5);
avg_t_pre_full=(mean(reshape(t_pre_full,trials,last-first+1)))';

h=loglog(N,t_ndft,'ko',avg_N,avg_t_ndft,'k-.',...
         N,t_horner,'k+',avg_N,avg_t_horner,'k--',...
         N,t_pre_full,'kx',avg_N,avg_t_pre_full,'k',...
         N,10^-8*N.^2,'k--',...
         N,2*10^-7*N.^2,'k-.');
set(h,'LineWidth',1.8); set(h,'MarkerSize',6); 
set(gca,'YTick',[10^-6,10^-4,10^-2,1,10^2]);
set(gca,'XTick',[10^2,10^3,10^4,10^5,10^6]);
set(gca,'FontSize',20);
axis([N(1),N(end),10^-7,10^2]);

print temp.eps -deps
!ps2pdf temp.eps ndft_fast1.pdf 
!rm temp.eps
