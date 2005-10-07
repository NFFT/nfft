clear;
n=128;
p=1;
a = 1:1500;
b = -1500:1500;
[A,B] = meshgrid(a,b); 
S = A.^2 + B.^2;
Z = exp(-pi^2*A*n^2./(4*p^2*S));
ZZ = 2*exp(-A*(2*p-1)^2/4).*(1+1./p*(2*p-1)*A) + sqrt(pi./sqrt(S))/p.*exp(-pi^2*A*n^2./(4*p^2*S)).*(1 + 2*S*p^2./(n*pi^2*A));

p =-20:2:2;
v = 10.^p;
figure(1);
[C,h]=contour(A,B,ZZ,v);
clabel(C,h);
figure(2);
[C,h]=contour(A,B,Z,v);
clabel(C,h);