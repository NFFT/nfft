function F = pfmatrix(N,theta,phi)
N
  [x,y] = meshgrid((-N):N,theta);
  f1 = exp(-i*x.*y);
  [x2,y2] = meshgrid((-N):N,phi);
  f2 = exp(-i*x2.*y2);
  F = [];
  for n=1:length(theta)
    F = [F;kron(f2(n,:),f1(n,:))];   
  end
