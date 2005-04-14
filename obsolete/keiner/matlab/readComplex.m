function x = readComplex(filename,n)
  f = fopen(filename);
  x = zeros(n,1);
  for k = 1:n
    x(k) = fscanf(f,'%f\n',1) + i * fscanf(f,'%f\n',1);
  end
  fclose(f);