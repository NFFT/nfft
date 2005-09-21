%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples of condition numbers of the kernel matrix K_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=400;
freq=(-N/2):(N/2-1);
Mset=[2,11,31];
Xtry=[20,200,200];

% Fejer kernel
F_omega_hat=2/N*(1-abs(2*freq'+1)/N);
zeta2=pi^2/6;
nom2=2*zeta2;

q=zeros(max(Xtry),2);
condK=zeros(max(Xtry),2);
for l=1:length(Mset)
  for j=1:Xtry(l)
    M=Mset(l);
    x=sort(rand(M,1)-0.5);
    xp=[x(end)-1;x;x(1)+1];
    q(j,l)=min(xp(2:end)-xp(1:end-1));

    A=exp(-2*pi*i*x*freq);
    condK(j,l)=cond(A*diag(F_omega_hat)*A');
  end
end

qreg=linspace(1/N,max(1./Mset),10000);
estimate_Lambda=1+nom2./(N^2*qreg.^2);
estimate_lambda=1-nom2./(N^2*qreg.^2);
estimate_cond=estimate_Lambda./estimate_lambda;

figure(1);
h=loglog(q(:,1),condK(:,1),'kx',q(:,2),condK(:,2),'ko',q(:,3),condK(:,3),'k+',qreg,estimate_cond,'k');
axis([1/N,max(1./Mset),0.9,10]);
set(h,'LineWidth',1.8); set(gca,'FontSize',25);
print test_cond1.eps -deps

% B-Spline kernel g_4
B_tilde_omega_hat=cubic_spline(freq'/N);
B_omega_hat=B_tilde_omega_hat+[B_tilde_omega_hat(2:end);0];
B_omega_hat=B_omega_hat/sum(B_omega_hat);
zeta4=pi^4/90;
nom4=2*zeta4*(zeta4+1)*8^4/((2*pi)^4-2*zeta4);

q=zeros(max(Xtry),2);
condK=zeros(max(Xtry),2);
for l=1:length(Mset)
  for j=1:Xtry(l)
    M=Mset(l);
    x=sort(rand(M,1)-0.5);
    xp=[x(end)-1;x;x(1)+1];
    q(j,l)=min(xp(2:end)-xp(1:end-1));

    A=exp(-2*pi*i*x*freq);
    condK(j,l)=cond(A*diag(B_omega_hat)*A');
  end
end

qreg=linspace(1/N,max(1./Mset),10000);
estimate_Lambda=1+nom4./(N^4*qreg.^4);
estimate_lambda=1-nom4./(N^4*qreg.^4);
estimate_cond=estimate_Lambda./estimate_lambda;

figure(2);
h=loglog(q(:,1),condK(:,1),'kx',q(:,2),condK(:,2),'ko',q(:,3),condK(:,3),'k+',qreg,estimate_cond,'k');
axis([1/N,max(1./Mset),0.9,10]);
set(h,'LineWidth',1.8); set(gca,'FontSize',25);
print test_cond2.eps -deps
