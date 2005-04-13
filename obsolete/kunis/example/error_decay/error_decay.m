% examples in shown are q=1/N, q=2/N, q=4/N, q=6/N
% output in error_decay.test and error_decay_1.eps, ...

format short
format compact
N=1000;
M=100;
l_max=100;

!make

q=6/N;
[X,q_X]=my_rand(M,q);
q_X

f_tilde=rand(M,2);
input_data=[X,f_tilde];
save input_data.dat -ascii -double -tabs input_data

!error_decay 1000 100 0 > output_data.dat
load output_data.dat
error0=output_data/output_data(1);

!error_decay 1000 100 1 > output_data.dat
load output_data.dat
error1=output_data/output_data(1);

!error_decay 1000 100 2 > output_data.dat
load output_data.dat
error2=output_data/output_data(1);

% estimate convergence rate numerically
n=find(error0<10^-10); n=n(1);
gamma_tilde_num0=(error0(n)/2)^(1/(n-1))
n=find(error1<10^-10); n=n(1);
gamma_tilde_num1=(error1(n)/2)^(1/(n-1))
n=find(error2<10^-10); n=n(1);
gamma_tilde_num2=(error2(n)/2)^(1/(n-1))

% theoretical convergence rate, see table 2.2
gamma0=(1+log(M/2))/(N*q_X);
if(gamma0<1)
  gamma_tilde0=(sqrt(1+gamma0)-sqrt(1-gamma0))/(sqrt(1+gamma0)+sqrt(1-gamma0))
else
  gamma_tilde0=inf
end;
gamma1=pi^2/(3*N^2*q_X^2);
if(gamma1<1)
  gamma_tilde1=(sqrt(1+gamma1)-sqrt(1-gamma1))/(sqrt(1+gamma1)+sqrt(1-gamma1))
else
  gamma_tilde1=inf
end;
gamma2=(2*2^4*pi^4)/(90*N^4*q_X^4);
if(gamma2<1)
  gamma_tilde2=(sqrt(1+gamma2)-sqrt(1-gamma2))/(sqrt(1+gamma2)+sqrt(1-gamma2))
else
  gamma_tilde2=inf
end;

l=1:20;
figure(1); 
h=semilogy(l-1,gamma_tilde_num0.^(l-1),'k',l-1,gamma_tilde0.^(l-1),'k--',l-1,error0(l),'ko');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',30);
%legend([texlabel('gamma='),sprintf('%1.1e',gamma_tilde_num0)],[texlabel('gamma='),sprintf('%1.1e',gamma_tilde0)]);
if(gamma_tilde0<inf)
  legend(sprintf('%1.1e',gamma_tilde_num0),sprintf('%1.1e',gamma_tilde0));
else
  legend(sprintf('%1.1e',gamma_tilde_num0));
end;
axis([l(1)-1,l(end)-1,10^-16,100]);
print error_decay_4_0.eps -deps

figure(2); 
h=semilogy(l-1,gamma_tilde_num1.^(l-1),'k',l-1,gamma_tilde1.^(l-1),'k--',l-1,error1(l),'kx');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',30);
if(gamma_tilde1<inf)
  legend(sprintf('%1.1e',gamma_tilde_num1),sprintf('%1.1e',gamma_tilde1));
else
  legend(sprintf('%1.1e',gamma_tilde_num1));
end;
axis([l(1)-1,l(end)-1,10^-16,100]);
print error_decay_4_1.eps -deps

figure(3); 
h=semilogy(l-1,gamma_tilde_num2.^(l-1),'k',l-1,gamma_tilde2.^(l-1),'k--',l-1,error2(l),'k+');
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',30);
if(gamma_tilde2<inf)
  legend(sprintf('%1.1e',gamma_tilde_num2),sprintf('%1.1e',gamma_tilde2));
else
  legend(sprintf('%1.1e',gamma_tilde_num2));
end;
axis([l(1)-1,l(end)-1,10^-16,100]);
print error_decay_4_2.eps -deps

