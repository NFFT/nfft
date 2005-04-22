function genData(filename,D,M,type,mode)
  
  f = fopen(filename,'w');

  fprintf(f,'%d\n',type);  
  fprintf(f,'%d\n',mode);
  
  fprintf(f,'%d\n',D);
  for i=1:D
    fprintf(f,'%17.16f\n%17.16f\n',rand-0.5,-0.5*rand);  
  end    
  
  fprintf(f,'%d\n',M);
  
  if (type == 0)
    last = (M+1)*(M+1);
  else
    last = D;
  end
  
  for k=1:last
    fprintf(f,'%17.16f\n',rand);
  end

  fclose(f);