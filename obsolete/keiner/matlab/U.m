function U = U(n,tau,L)
  C3 = DCT3matrix(2^(tau+1));
  ZP = Z(tau);
  D3 = kron(sparse([1,0;0,1;1,0;0,1]),C3)*ZP;
  x = cos((2*(0:2^(tau+1)-1)+1)*pi/(2^(tau+2)));
  c = (2^(tau+1))*L+1;
  k = 2^(tau)-1;
  if n == 0
    %k
    %n
    %c
  end
  %if n == 0
    [alpha,beta,gamma] = alcoeff(c,n);
    PV = sparse(diag([gamma*aleval(k-1,n,c+1,x),gamma*aleval(k,n,c+1,x),aleval(k,n,c,x),aleval(k+1,n,c,x)]));
  %else
  %  PV = diag([gamma(c,n)*P(k-1,n,c+1,x,0),gamma(c,n)*P(k,n,c+1,x,0),P(k,n,c,x,0),P(k+1,n,c,x,0)]);
  %end
  %P = diag(ones(2^(tau+3),1));
  S = kron(speye(2),[speye(2^(tau+1)),speye(2^(tau+1))]);
  C2 = DCT2matrix(2^(tau+1));
  D = speye(2^(tau+1));
  D(1,1) = 0.5;
  Dn = (2^(-tau))*speye(2^(tau+1));
  D2 = kron(speye(2),Dn*D*C2);
  size(D2);
  size(S);
  size(PV);
  size(D3);
  U = D2*S*PV*D3;
  if n == 0 
    %D2
    %S
    %diag(PV)
    %D3
  end
