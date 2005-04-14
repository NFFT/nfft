function alpha = alpha(k,n)
  if k < n
    alpha = (-1)^(k+1);
  else
    alpha = (2*k+1)/(sqrt((k-n+1)*(k+n+1)));
  end