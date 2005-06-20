function [M] = construct_readout_time( M, Te, arms, increment )

out=zeros(M,1);

for j=1:M,
  out(j) =   mod(j-1,M/arms)*increment + Te;
end

%0.00402542373

save readout_time.dat -ascii out
