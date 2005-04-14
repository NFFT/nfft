function gamma = gamma(k,n)
if k <= n
  gamma = 0;
else
  gamma = -(sqrt((k-n)*(k+n)))/(sqrt((k-n+1)*(k+n+1)));
end