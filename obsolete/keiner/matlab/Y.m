function y = Y(n,theta,phi)
  p = legendre(abs(n),cos(theta),'sch'); 
  if n>0
    p(2:end,:) = (1/sqrt(2))*p(2:end,:);
  end
  [A,B] = meshgrid(phi,(-n:1:n));
  y = [flipud(p(2:end,:));p] .* exp(i*A.*B);
 