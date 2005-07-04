function generateNodes(d_phi,d_theta,alpha)
  D = round(d_theta*d_phi)
  P = randperm(D);
  D2 = round(alpha*D);
  file = fopen('nodes.dat','w');
  fprintf(file,'%d\n%d\n%d\n',d_theta,d_phi,D2);
  fprintf(file,'%d\n',P(1:D2)-1);
  fclose(file);
