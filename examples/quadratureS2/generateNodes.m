function generateNodes(M)
% ACCURACY Generate Gauss-Legendre weights.

prefix = 'gl';
suffix = '.dat';

for k = M
  % Generate Gauss-Legendre nodes in latitudinal direction.
  [theta,w] = lgwt(k+1,-1,1);
  theta = (1/(2*pi))*acos(theta);
  
  % Generate nodes in colatitudinal direction.
  phi = (1/(2*pi))*(0:2*k+1).*(pi/(k+1));
  
  % Open file for writing. */
  file = fopen([prefix,int2str(k),suffix],'w');
  
  % Write data to file.
  fprintf(file,'%.30f\n',theta);
  fprintf(file,'%.30f\n',phi);
  fprintf(file,'%.30f\n',w);
  
  % Close datafile.
  fclose(file);
end
