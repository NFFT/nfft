function [Y,F,a,theta,phi,N,t] = readData(filename)
  f = fopen(filename,'r');
  
  type = fscanf(f,'%d\n',1);
  mode = fscanf(f,'%d\n',1);
  
  D = fscanf(f,'%d\n',1);
  phi = zeros(D,1);
  theta = zeros(D,1);
  for i=1:D
    phi(i) = fscanf(f,'%f\n',1);
    theta(i) = fscanf(f,'%f\n',1);    
  end    
  
  M = fscanf(f,'%d\n',1);
  t = ceil(log2(M));
  N = 2^t;
  Y = sfmatrix(M,2*pi*theta',2*pi*phi','semi');
  F = pfmatrix(N,2*pi*theta',2*pi*phi');
  
  if (type == 0)
    last = (M+1)^2;
  else
    last = D;
  end
  
  a = zeros(last,1);
  for i=1:last
    a(i) = fscanf(f,'%f\n',1);
  end
  
  fclose(f);
