function [Y,F] = Ymatrix(M,N,theta,phi)
  Y = [];
  F = [];
  for n=0:M
    Y = [Y;Yf(n,theta,phi)];
  end
  Y = Y';
 
  [x,y] = meshgrid((-N):N,theta);
  f1 = exp(-i*x.*y);
  [x2,y2] = meshgrid((-N):N,phi);
  f2 = exp(-i*x2.*y2);
  for n=1:length(theta)
    F = [F;kron(f2(n,:),f1(n,:))];   
  end
  %exp(-2*pi*theta)
  