function g=sobolev_weight(z,alpha,beta,gamma)

if(beta==0)
  g=(gamma+abs(z).^(2*alpha)).^(-1);
else
  g=(gamma+abs(z).^(2*alpha)).^(-1) .* (1/4-z.^2).^beta * 4^beta;
end
