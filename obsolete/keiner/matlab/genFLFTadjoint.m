function T = genFLFT(M,t)
N = 2^t;
for n=-M:M
  nleg = abs(n); 
  % Generate T_0
  e = zeros(2*N,1);
  
  [alpha,beta,gamma] = alcoeff(N-1,nleg);
  e(2*N) = alpha;
  e(2*N-1) = beta;
  e(2*N-3) = gamma;
  
  %fprintf('n = %d\n',n);
  %fprintf('alpha = %f\n',alpha(N-1,nleg));
  %fprintf('beta = %f\n',beta(N-1,nleg));
  %fprintf('gamma = %f\n',gamma(N-1,nleg));
  
  T0 = [kron(eye(N),[1;0]),e];
  T{n+M+1,1} = T0;
  % Generate T_Tau
  for tau = 1:t-1
    Ttau = zeros(2*N);
    for L=0:(2^(t-tau-1)-1)
      %fprintf('U(nleg = %d, tau = %d, L = %d)\n',nleg,tau,L);  
      Ttau(L*(2^(tau+2))+1:(L+1)*(2^(tau+2)),L*(2^(tau+2))+1:(L+1)*(2^(tau+2))) = [Z(tau),U(nleg,tau,L)];
    end
    T{n+M+1,tau+1} = Ttau;
  end
  % Generate T_t
  nleg = abs(n);  
  gamman = sqrt(factorial(2*nleg))/(2^nleg*factorial(nleg));
  if mod(n,2) == 1
    Tt = [eye(N+1),eye(N+1)];
  else
    W = 0.5*(diag(ones(1,N),1) + diag(ones(1,N),-1));
    W (2,1) = 1;
    if n == 0
      Tt = [eye(N+1),W];
    else
      Tt = [eye(N+1),eye(N+1)-W];
    end
  end
  T{n+M+1,t+1} = (-1)^n * gamman*Tt(1:N+1,[1:N,N+2:end-1]);
end  
