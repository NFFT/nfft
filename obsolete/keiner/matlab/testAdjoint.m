function err = testAdjoint(D,M,mode)
% mode = {0,1,2,3,4,5}
% 0 = Compare forward implementation with bare Y
% 1 = Compare forward implementation with factorized Y
% 2 = Compare bare Y with factorized Y
% 3 = Compare adjoint implementation with bare Y^H
% 4 = Compare adjoint implementation with factorized Y^H
% 5 = Compare bare Y^H with factorized Y^H

% Some constants.
filename = 'temp';
insuffix = '.in';
outsuffix = '.out';

% Generate new data set.
if (mode == 0 || mode == 1 || mode == 2)
  genData([filename,insuffix],D,M,0);
else
  genData([filename,insuffix],D,M,1);
end

% Generate Y, the Fourier matrix F and the data vector a.
[Y,F,a,theta,phi,N,t] = readData([filename,insuffix]);

% Compute decomposed Y, if neccesary.
if (mode == 1 || mode == 2 || mode == 4 || mode == 5)    
  A = F';
  A = [sparse((2*M+1)*(2*N+1),(N-M)*(2*N+1)),speye((2*M+1)*(2*N+1)),sparse((2*M+1)*(2*N+1),(N-M)*(2*N+1))] * A;

  % Adjoint Cheb2Exp
  S = 0.5*i*(sparse(diag(ones(2*N,1),1)) - sparse(diag(ones(2*N,1),-1)));
  X = speye((2*M+1)*(2*N+1));
  for n=-M:M
    if mod(n,2) == 1
      X((n+M)*(2*N+1)+1:(n+M+1)*(2*N+1),(n+M)*(2*N+1)+1:(n+M+1)*(2*N+1)) = S;  
    end
  end
  A = X*A;
  
  C = 0.5*speye(N);
  T = sparse([zeros(1,N),1,zeros(1,N);fliplr(C),zeros(N,1),C]);
  W = sparse(kron(eye(2*M+1),T));
  A = W*A;

  % Adjoint FLFT
  T = genFLFTadjoint(M,t);
  V = sparse((2*M+1)*(N+1),(2*M+1)*(N+1));
  for n=-M:M
    Ttau = T{n+M+1,t+1}';
    for j=t-1:-1:1
      size(T{n+M+1,j+1}');
      Ttau = T{n+M+1,j+1}' * Ttau;
    end  
    Ttau = T{n+M+1,1}' * Ttau;
    V((n+M)*(N+1)+1:(n+M+1)*(N+1),(n+M)*(N+1)+1:(n+M+1)*(N+1)) = Ttau;
  end

  A = V*A;

  Z = sparse((M+1)^2,(2*M+1)*(N+1));
  row = 1;
  for k=0:M 
    for n = -k:k
     Z(row,(n+M)*(N+1)+1+k) = 1; 
     row = row + 1;
    end
  end

  A = Z*A;  
end

if (mode == 0 || mode == 1)
!../examples/nfsft < temp.in > temp.out
% Read forward transformed data.
ri = readComplex([filename,outsuffix],size(Y,1));  
else
  if  (mode == 3 || mode == 4)
    % Read adjoint transformed data.
    ri = readComplex([filename,outsuffix],size(Y',1));
  end
end    

if (mode == 0 || mode == 2)
% Bare forward transform.
rb = Y*a;  
else
  if  (mode == 3 || mode == 5)
    % Bare adjoint transform.
    rb = Y'*a;
  end
end    

if (mode == 1 || mode == 2)
% Decomposed forward transform.
rd = A'*a;  
else
  if  (mode == 4 || mode == 5)
    % Decomposed adjoint transform.
    rd = A*a;
  end
end    

switch mode
  case 0
    fprintf('Forward implementation vs. bare Y:\n');
    norm(ri-rb,1)/norm(rb,1)
    %image(abs(ri-rb));
    %colorbar;
  case 1
    fprintf('Forward implementation vs. decomposed Y:\n');
    norm(ri-rd,1)/norm(rd,1)
    %image(abs(ri-rd));
    %colorbar;
  case 2
    fprintf('Bare Y vs. decomposed Y:\n');
    norm(rd-rb,1)/norm(rb,1)
    %image(abs(rd-rb));
    %colorbar;
  case 3
    fprintf('Adjoint implementation vs. bare Y^H:\n');
    norm(ri-rb,1)/norm(rb,1)
  case 4
    fprintf('Adjoint implementation vs. decomposed Y^H:\n');
    norm(ri-rd,1)/norm(rd,1)
  case 5
    fprintf('Bare Y^H vs. decomposed Y^H:\n');
    norm(rd-rb,1)/norm(rb,1)
end

return;

Y = Y';
Y(1:10,1:5)
A(1:10,1:5)
norm(A-Y,1)/norm(Y,1)

return;

r = readComplex([filename,outsuffix],length(z));
d = z-r;
norm(d,1)



r2 = Y'*a;
[r,r2];
g = r-r2;
d = zeros(M+1,2*M+1);
ll = 1;
for k=0:M
  d(k+1,M+1-k:M+1+k) = g(ll:ll+2*k)';
  ll = ll + 2*k+1;
end
image(abs(d));
err = norm(d,1)/norm(r2);
return;

