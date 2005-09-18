function p=my_legendre(k,x)

pp=legendre(k,x);
if(k==0)
  p=pp;
else
  p=squeeze(pp(1,:,:));
end;
