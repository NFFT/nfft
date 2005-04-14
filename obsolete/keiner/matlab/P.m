function y = P(k,n,c,x,debug)
if debug == 1
  fprintf('k = %d\n',k);
end

  if k == -1
    y = zeros(size(x));
  else
    if k == 0
      y = ones(size(x));
    else
      a = ones(size(x));
      b = zeros(size(x));
      for j = k:-1:2
        a_old = a;
        if debug == 1
          fprintf('j = %d\n',j-1+c);  
          alpha(j-1+c,n)
          beta(j-1+c,n)
          gamma(j-1+c,n)
        end  
        a = b + alpha(j-1+c,n)*x.*a_old + beta(j-1+c,n)*a_old;
        b = gamma(j-1+c,n)*a_old;
      end    
      if debug == 1
        alpha(c,n)
        beta(c,n)
        end  
      y = alpha(c,n)*x.*a + beta(c,n)*a + b;
    end
  end