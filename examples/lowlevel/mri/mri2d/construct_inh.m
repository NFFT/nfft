function [] = construct_inh( N )
out=zeros(N*N,1);

for k=1:N,
  for j=1:N,
    out((k-1)*N+j) = (((k-N/2)/(N/sqrt(2)))^2+((j-N/2)/(N/sqrt(2)))^2)/2 -0.5/2;
  end
end

save inh.dat -ascii out
