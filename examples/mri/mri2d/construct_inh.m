function [] = construct_inh( N )
out=zeros(N*N,1);

for k=1:N,
  for j=1:N,
    out((k-1)*N+j) = sqrt((((k-N/2)/(N/sqrt(2)))^2+((j-N/2)/(N/sqrt(2)))^2)) -0.5;
  end
end

imagesc(reshape(out,N,N));

save inh.dat -ascii out
