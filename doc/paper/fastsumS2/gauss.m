sigma=0.1;
n=20;

g=zeros(n,1);

g(0+1)=4*pi*sigma^(-1)*exp(-2*sigma)*sinh(sigma)*cosh(sigma);
g(1+1)=pi*sigma^(-2)*exp(-2*sigma)*(2*sigma*cosh(2*sigma)-sinh(2*sigma));

for k=1:n-2
 g(k+1+1)=g(k-1+1) - (2*k+1)/(2*sigma) * g(k+1);
 if g(k+1+1)<0;
  g=g(1:k+1);
  break;
 end;
end;

log(g)
