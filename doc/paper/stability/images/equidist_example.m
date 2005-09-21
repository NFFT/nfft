M=100;
N_min=M;
N_max=6*M;

cond_D=zeros(N_max,1);
cond_F=zeros(N_max,1);

for N=N_min:2:N_max
  freq=-N/2:(N/2-1);

  w_hat_D=ones(size(freq))/N;	
  w_hat_F=2/N*(1-abs(2*freq+1)/N);

  lambda_D=zeros(M,1);
  lambda_F=zeros(M,1);	
  for j=1:M
    lambda_D(j)=M*sum(w_hat_D(j:M:end));
    lambda_F(j)=M*sum(w_hat_F(j:M:end));
  end;	

  cond_D(N)=max(lambda_D)/min(lambda_D);
  cond_F(N)=max(lambda_F)/min(lambda_F);	 
end

figure(1);
N=N_min:2:N_max;
Bound_cond_F=(N.^2+M.^2)./(N.^2-M.^2);
h=loglog(N,cond_D(N),'k-.',N,cond_F(N),'k',N,Bound_cond_F,'k--');
axis([N_min,N_max,1,sqrt(M)]);
set(h,'LineWidth',1.8); set(gca,'FontSize',25);
set(gca,'XTick',[N_min,2*N_min,N_max]);
set(gca,'YTick',[1,2,10]);

print equidist_example1.eps -deps
