f = fopen('convolution.dat','r');
D = fscanf(f,'%d\n',1);
D_PHI = fscanf(f,'%d\n',1);
D_THETA = fscanf(f,'%d\n',1);
x = zeros(D_THETA,D_PHI);
y = zeros(D_THETA,D_PHI);
z = zeros(D_THETA,D_PHI);
for d_theta = 1:D_THETA
  for d_phi = 1:D_PHI
    phi = fscanf(f,'%E ',1);
    phi = 2*pi*phi;
    theta = fscanf(f,'%E ',1);
    theta = 2*pi*theta;
    r = 1 + fscanf(f,'%E\n',1);
    x(d_theta,d_phi) = r*sin(theta)*cos(phi);
    y(d_theta,d_phi) = r*sin(theta)*sin(phi);
    z(d_theta,d_phi) = r*cos(theta);
  end
end
fclose(f);
surf(x,y,z);
