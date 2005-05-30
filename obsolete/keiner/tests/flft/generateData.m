function generateNodes(M)

prefix = 'gl';
suffix = '.dat';

for k = 1:size(M,1)
  % Open file for writing
  file = fopen([prefix,'_',int2str(M(k,1)),'_',int2str(M(k,2)),suffix],'w');
  
  % Write data to file.
  fprintf(file,'%d\n',M(k,1));
  fprintf(file,'%d\n',M(k,2));
  
  for i = M(k,2):M(k,1)
    fprintf(file,'%.30f\n%.30f\n',rand(1),rand(1));
  end
  
  % Close datafile.
  fclose(file);
end
