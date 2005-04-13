load timing_1d.dat
timing_1d(timing_1d<0.01)=inf;
N=timing_1d(:,1);
N=N(1:end/3);
anz=size(N,1);
fast1=timing_1d(1:anz,2);
fast2=timing_1d(anz+1:2*anz,2);
fast3=timing_1d(2*anz+1:3*anz,2);
slow=timing_1d(1:anz,3);

h=loglog(N,slow,'k',N,slow,'ko',...
	 N,fast1,'k',N,fast1,'kx',...
	 N,fast2,'k',N,fast2,'k+',...
	 N,fast3,'k',N,fast3,'k^');
axis([2^10 2^20 10^-1 10^2]);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22); set(gca,'XTick',N(1:2:end)); print timing1.eps -deps


load timing_2d.dat
timing_2d(timing_2d<0.01)=inf;
N=timing_2d(:,1);
fast=timing_2d(:,2);
slow=timing_2d(:,3);

h=loglog(N,slow,'k',N,slow,'ko',...
	 N,fast,'k',N,fast,'k+');
axis([2^4 2^10 10^-1 10^2]);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22); set(gca,'XTick',N(1:2:end)); print timing2.eps -deps

load timing_3d.dat
timing_3d(timing_3d<0.01)=inf;
N=timing_3d(:,1);
fast=timing_3d(:,2);
slow=timing_3d(:,3);

h=loglog(N,slow,'k',N,slow,'ko',...
	 N,fast,'k',N,fast,'k+');
axis([2^4 2^6 10^-1 10^2]);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22); set(gca,'XTick',N(1:2:end)); print timing3.eps -deps


load timing_2d_fast1.dat
load timing_2d_fast2.dat
load timing_2d_fast3.dat
timing_2d_fast1(timing_2d_fast1<0.01)=inf;
timing_2d_fast2(timing_2d_fast2<0.01)=inf;
timing_2d_fast3(timing_2d_fast3<0.01)=inf;
N=timing_2d_fast1(:,1);
fast1=timing_2d_fast1(:,2:4);
fast2=timing_2d_fast2(:,2:4);
fast3=[timing_2d_fast3(:,2:4);inf*ones(2,3)];

h=loglog(N,fast1(:,1),'k',N,fast1(:,1),'kx',...
	 N,fast2(:,1),'k',N,fast2(:,1),'k+',...
	 N,fast3(:,1),'k',N,fast3(:,1),'k^');
axis([2^7 2^10 10^-2 10^0]);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22); set(gca,'XTick',N(1:2:end)); print timing4.eps -deps

h=loglog(N,fast1(:,2),'k',N,fast1(:,2),'kx',...
	 N,fast2(:,2),'k',N,fast2(:,2),'k+',...
	 N,fast3(:,2),'k',N,fast3(:,2),'k^');
axis([2^7 2^10 10^-2 10^1]);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22); set(gca,'XTick',N(1:2:end)); print timing5.eps -deps

h=loglog(N,fast1(:,3),'k',N,fast1(:,3),'kx',...
	 N,fast2(:,3),'k',N,fast2(:,3),'k+',...
	 N,fast3(:,3),'k',N,fast3(:,3),'k^');
axis([2^7 2^10 10^-1 10^2]);
set(h,'LineWidth',2.0); set(h,'Markersize',10);
set(gca,'FontSize',22); set(gca,'XTick',N(1:2:end)); print timing6.eps -deps