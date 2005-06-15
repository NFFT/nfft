function generateNodes(n)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes her
[theta,phi] = heal_coord(n);
file = fopen('healpix.dat','w');
fprintf(file,'%d\n',n);
for i = 1:length(theta)
  theta(i) = (0.5/pi)*theta(i);
  phi(i) = (0.5/pi)*phi(i);
  if (phi(i) >= 0.5)
    phi(i) = phi(i) - 1.0;
  end
  fprintf(file,'%.16f %.16f\n',theta(i),phi(i));
end
fclose(file);
