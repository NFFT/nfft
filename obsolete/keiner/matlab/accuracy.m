function err = accuracy(M,rep)
% ACCURACY Test accuracy of explicit discrete spherical Fourier transform.

err = zeros(size(M));

i = 1;
for j = M
  % Generate Gauss-Legendre nodes in latitudinal direction.
  [theta,w] = lgwt(j+1,-1,1);
  theta = acos(theta);

  % Generate nodes in colatitudinal direction.
  phi = (0:2*j+1).*(pi/(j+1));

  w = repmat(w,1,length(phi));
  w = w(:);
  w = sparse(1:length(w),1:length(w),w);

  [phi,theta] = meshgrid(phi,theta);

  theta = theta(:)';
  phi = phi(:)';

  Y = sfmatrix(j,theta,phi,'norm');
  for d = 1:rep
    a = rand((j+1)*(j+1),1);
    z = (1/(2*j+2))*Y'*w*Y*a;
    err(i) = err(i) + max(abs(z-a));
  end
  err(i) = err(i)/rep;
  i = i + 1;
  fprintf('%d\n',j);
end