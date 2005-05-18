function [err,r1,r2] = testNFSFT(D,M,test,plot)
% TEST
% [ERR,R1,R2] = TEST(D,M,TEST)
% TEST \in {0,1,2,3,4,5,6,7,8,9}
% 0 = Compare fast forward implementation with bare Y
% 1 = Compare direct forward implementation with bare Y
% 2 = Compare fast forward implementation with factorized Y
% 3 = Compare direct forward implementation with factorized Y
% 4 = Compare bare Y with factorized Y
% 5 = Compare fast adjoint implementation with bare Y^H
% 6 = Compare direct adjoint implementation with bare Y^H
% 7 = Compare fast adjoint implementation with factorized Y^H
% 8 = Compare direct adjoint implementation with factorized Y^H
% 9 = Compare bare Y^H with factorized Y^H

% Argument checks
if (isnumeric(D) == 0 || isequal(size(D),[1,1]) == false || D <= 0)
  error('testNFSFT: Parameter D must be a positive natural number!');
end

if (isnumeric(M) == 0 || isequal(size(M),[1,1]) == false || M < 0)
  error('testNFSFT: Parameter M must be a natural number!');
end

if (isnumeric(test) == 0 || isequal(size(test),[1,1]) == false || test < 0 || test > 9)
  error('testNFSFT: Parameter TEST must be a natural number from [0,9]!');
end

% Some constants.
filename = 'temp';
insuffix = '.in';
outsuffix = '.out';

% Reset random number generator.
rand('state',0);

if (test >= 0 && test <= 4)
  type = 0
else
  type = 1
end

if ((test >= 0 && test <= 3) || (test >= 5 && test <= 8))
  usec = true
else
  usec = false
end

if (sum(test == [0,2,4,5,7]) > 0)
  mode = 0;
end
	
if (sum(test == [1,3,6,8,9]) > 0)
  mode = 1;
end

% Generate new data set.
genData([filename,insuffix],D,M,type,mode);

% Generate Y, the Fourier matrix F and the data vector a.
[Y,F,a,theta,phi,N,t] = readData([filename,insuffix]);

% Compute decomposed Y, if neccesary.
if (sum(test == [2,3,4,7,8,9]) > 0)    
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
  t
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

if (usec == true)
  !../examples/wrapper/wrapper < temp.in > temp.out
	if (type == 0)
	  % Read forward transformed data.
    ri = readComplex([filename,outsuffix],size(Y,1));  
  else
    % Read adjoint transformed data.
    ri = readComplex([filename,outsuffix],(M+1)^2);
  end
end    
	
if (sum(test == [0,1,4]) > 0)
% Bare forward transform.
rb = Y*a; 
else
  if (sum(test == [5,6,9]) > 0)
    % Bare adjoint transform.
    rb = Y'*a;
  end
end    

if (sum(test == [2,3,4]) > 0)
% Decomposed forward transform.
rd = A'*a;  
else
  if  (sum(test == [7,8,9]) > 0)
    % Decomposed adjoint transform.
    rd = A*a;
  end
end    

switch test
  case 0
    fprintf('Fast forward implementation vs. bare Y:\n');
	  r1 = rb;
	  r2 = ri;
  case 1
    fprintf('Direct forward implementation vs. bare Y:\n');
	  r1 = rb;
	  r2 = ri;
  case 2
    fprintf('Fast forward implementation vs. decomposed Y:\n');
	  r1 = rd;
	  r2 = ri;
  case 3
    fprintf('Slow forward implementation vs. decomposed Y:\n');
	  r1 = rd;
	  r2 = ri;
  case 4
    fprintf('Bare Y vs. decomposed Y:\n');
	  r1 = rb;
	  r2 = rd;
  case 5
    fprintf('Fast adjoint implementation vs. bare Y^H:\n');
	  r1 = rb;
	  r2 = ri;
  case 6
    fprintf('Direct adjoint implementation vs. bare Y^H:\n');
	  r1 = rb;
	  r2 = ri;
  case 7
    fprintf('Fast adjoint implementation vs. decomposed Y^H:\n');
	  r1 = rd;
	  r2 = ri;
  case 8
    fprintf('Direct adjoint implementation vs. decomposed Y^H:\n');
	  r1 = rd;
	  r2 = ri;
  case 9
    fprintf('Bare Y^H vs. decomposed Y^H:\n');
	  r1 = rb;
	  r2 = rd;
end

err = norm(r1-r2,1)/norm(r1,1);

fprintf('|r1|: max = %f, min = %f, NaN = %d\n', max(abs(r1)), min(abs(r1)), sum(sum(isnan(r1))));
fprintf('|r2|: max = %f, min = %f, NaN = %d\n', max(abs(r2)), min(abs(r2)), sum(sum(isnan(r2))));
fprintf('||r1-r2||_1/||r1||_1 = %17.16f\n', err);

if (plot == 1)
  if (test <= 4)
    rm = abs(r1-r2);
    figure;
    spy(rm>=1E-7);
    colorbar;
    figure
    imagesc(rm);
    colorbar;
  else
    r = abs(r1-r2);
    rm = zeros(M+1,2*M+1);
    index = 1;
    for k = 0:M
      for n = -k:k
      rm(k+1,n+M+1) = r(index);
      index = index + 1;
      end
    end
    figure;
    spy(rm>=1E-7);
    colorbar;
    figure
    imagesc(rm);
    colorbar;
  end
end
