%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples of condition numbers of the kernel matrix K_N for jitter error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberAvg=100;
Mmax=100;
Mset=1:Mmax;
releps=0.1;

zeta2=pi^2/6;

condK=zeros(Mmax,1);
condK2=zeros(Mmax,1);
estimate_cond=zeros(Mmax,1);
for a=1:numberAvg
 for M=Mset
  xeq=(-0.5:(1/M):(0.5-1/M))';
  x=xeq+releps*1/M*rand(M,1);
  N=2*ceil(6*M/2);
  freq=(-N/2):(N/2-1);

  % Fejer kernel
  F_omega_hat=2/N*(1-abs(2*freq'+1)/N);
 
  A=exp(-2*pi*i*x*freq);
  %condK(M)=condK(M)+cond(A*diag(F_omega_hat)*A');
  condK(M)=max([condK(M),cond(A*diag(F_omega_hat)*A')]);
  %condK2(M)=condK2(M)+cond(A*A');
  condK2(M)=max([condK2(M),cond(A*A')]);

  estimate_Lambda=1+2*zeta2*M^2 / (N^2 * (1-releps)^2);
  estimate_lambda=1-2*zeta2*M^2 / (N^2 * (1-releps)^2);
  estimate_cond(M)=estimate_Lambda/estimate_lambda;
 end
end

%condK=condK/numberAvg;
%condK2=condK2/numberAvg;

estimate_cond(estimate_cond<1)=nan;
figure(1);
h=plot(Mset,condK(Mset),'kx',Mset,condK2(Mset),'k+',Mset,estimate_cond(Mset),'k--');
%axis([0,1/M,0.9,xx]);
set(h,'LineWidth',1.8); set(h,'Markersize',8); set(gca,'FontSize',25);
print test_cond1.eps -deps


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples of condition numbers of the kernel matrix K_N for jitter error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=100;
freq=(-N/2):(N/2-1);
M=20;
xeq=(-0.5:(1/M):(0.5-1/M))';
neps=100;
maxeps=linspace(0,1/M,neps);

% Fejer kernel
F_omega_hat=2/N*(1-abs(2*freq'+1)/N);
zeta2=pi^2/6;

condK=zeros(neps,1);
for j=1:neps
  x=xeq+maxeps(j)*rand(M,1);
  A=exp(-2*pi*i*x*freq);
  condK(j)=cond(A*diag(F_omega_hat)*A');
end

estimate_Lambda=1+2*zeta2/ N^2 * (1/M-maxeps).^(-2);
estimate_lambda=1-2*zeta2/ N^2 * (1/M-maxeps).^(-2);
estimate_lambda(estimate_lambda<0)=0;
estimate_cond=estimate_Lambda./estimate_lambda;

figure(1);
h=loglog(maxeps,condK,maxeps,estimate_cond);
axis([0,1/M,0.9,100]);
%set(h,'LineWidth',1.8); set(gca,'FontSize',25);
print test_cond1.eps -deps
