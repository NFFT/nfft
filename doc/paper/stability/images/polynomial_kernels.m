%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples of localised trigonometric polynomial kernels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=20;
freq=(-N/2):(N/2-1);

% Dirichlet kernel
D_omega_hat=ones(N,1)/N;

% Fejer kernel
F_omega_hat=2/N*(1-abs(2*freq'+1)/N);

% B-Spline kernel g_4
B_tilde_omega_hat=cubic_spline(freq'/N);
B_omega_hat=B_tilde_omega_hat+[B_tilde_omega_hat(2:end);0];
B_omega_hat=B_omega_hat/sum(B_omega_hat);

alpha=1;
beta=2;
gamma=10^-4;
S_tilde_omega_hat=sobolev_weight(freq'/N,alpha,beta,gamma);
S_omega_hat=S_tilde_omega_hat+[S_tilde_omega_hat(2:end);0];
S_omega_hat=S_omega_hat/sum(S_omega_hat);

x=linspace(-0.5,0.5,1000)';
D_0=exp(-2*i*pi*x*freq)*D_omega_hat;
F_0=exp(-2*i*pi*x*freq)*F_omega_hat;
B_0=exp(-2*i*pi*x*freq)*B_omega_hat;
S_0=exp(-2*i*pi*x*freq)*S_omega_hat;

figure(1); h=plot(x,real(D_0),'k'); axis([-0.5,0.5,-0.3,1.1]); 
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels1.eps -deps

figure(2); h=plot(x,real(F_0),'k'); axis([-0.5,0.5,-0.3,1.1]);
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels2.eps -deps

figure(3); h=plot(x,real(B_0),'k'); axis([-0.5,0.5,-0.3,1.1]);
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels3.eps -deps

figure(4); h=plot(x,real(S_0),'k'); axis([-0.5,0.5,-0.3,1.1]); 
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels4.eps -deps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decay of the B-Spline kernel B_{4,N} and the Jackson kernel J_{4,N}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma=21;

N=4*(sigma-1);

freq=(-N/2):(N/2-1);

B_tilde_omega_hat=cubic_spline(freq'/N);
B_omega_hat=B_tilde_omega_hat+[B_tilde_omega_hat(2:end);0];
B_omega_hat=B_omega_hat/sum(B_omega_hat);

x=linspace(0,0.5,1000)';
B_0=exp(-2*i*pi*x*freq)*B_omega_hat;
J_0=(1+exp(-2*pi*i*x))/(2*sigma^4).*(sin(sigma*pi*x)./sin(pi*x)).^4;

figure(5);
h=semilogy(x,abs(B_0),'k', x,16/3*abs(N*x).^(-4),'k--');
%h=plot(x,abs(B_0)./(16/3*abs(N*x).^(-4)),'k--');
axis([0,0.5,10^-8,1]);
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels5.eps -deps

sigma=20;
N=4*(sigma-1)+2;

freq=(-N/2):(N/2-1);

%B_tilde_omega_hat=cubic_spline(freq'/N);
%B_omega_hat=B_tilde_omega_hat+[B_tilde_omega_hat(2:end);0];
%B_omega_hat=B_omega_hat/sum(B_omega_hat);

x=linspace(0,0.5,1000)';
%B_0=exp(-2*i*pi*x*freq)*B_omega_hat;
J_0=(1+exp(-2*pi*i*x))/(2*sigma^4).*(sin(sigma*pi*x)./sin(pi*x)).^4;

figure(6);
h=semilogy(x,abs(J_0),'k',x,16*abs(N*x).^-4,'k--',x,abs(B_0),'k--');
%				h=plot(x,abs(J_0./B_0));
axis([0,0.5,10^-8,1]);
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels6.eps -deps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sobolev kernel and its weight function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=40;
freq=(-N/2):(N/2-1);

alpha=1;
beta=2;
gamma=10^-4;
S_tilde_omega_hat=sobolev_weight(freq'/N,alpha,beta,gamma);
S_omega_hat=S_tilde_omega_hat+[S_tilde_omega_hat(2:end);0];
S_omega_hat=S_omega_hat/sum(S_omega_hat);

S_tilde_omega_hat2=sobolev_weight(freq'/N,alpha,0,gamma);
S_omega_hat2=S_tilde_omega_hat2+[S_tilde_omega_hat2(2:end);0];
S_omega_hat2=S_omega_hat2/sum(S_omega_hat2);

x=linspace(0,0.5,1000)';
S_0=exp(-2*i*pi*x*freq)*S_omega_hat;
S_02=exp(-2*i*pi*x*freq)*S_omega_hat2;

figure(7); h=plot(x,real(S_0),'k',x,real(S_02),'k--'); axis([0,0.5,-0.3,1.1]); 
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels7.eps -deps

z=linspace(0,0.5,1000)';
figure(8); h=semilogy(z,sobolev_weight(z,alpha,beta,gamma),'k',z,z.^(-2*alpha),'k--'); axis([0,0.5,10^-4,10^4]); 
set(h,'LineWidth',2.0); set(gca,'FontSize',30);
print polynomial_kernels8.eps -deps
