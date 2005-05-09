function genData(filename,D,M,type,mode)
% GENDATA
%   GENDATA(FILENAME,D,M,TYPE,MODE)
%
  
  f = fopen(filename,'w');

  fprintf(f,'%d\n',type);  
  fprintf(f,'%d\n',mode);
  
  fprintf(f,'%d\n',D);
  for i=1:D
    fprintf(f,'%17.16f\n%17.16f\n',rand-0.5,0.5*rand);  
    %fprintf(f,'%17.16f\n%17.16f\n',0,0.25);  
  end    
  
  fprintf(f,'%d\n',M);
  
  if (type == 0)
    last = (M+1)*(M+1);
  else
    last = D;
  end
  
  skip1 = 0;
  skip2 = last;
  for k = 1:skip1
    fprintf(f,'%17.16f\n',0.0);
  end
  for k = skip1+1:skip2
    fprintf(f,'%17.16f\n',rand);
    %fprintf(f,'%17.16f\n',1);
  end
  for k = skip2+1:last
    fprintf(f,'%17.16f\n',0.0);
  end

  fclose(f);
