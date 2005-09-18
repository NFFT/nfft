maxM=19;
L=400;
D=L;

h=0.1;

X=2*rand(D,3)-1;
X=X./(sqrt(sum(X.^2,2))*ones(1,3));
Y=2*rand(L,3)-1;
Y=Y./(sqrt(sum(Y.^2,2))*ones(1,3));

Q=(1-h^2)./(4*pi*(1-2*h*X*Y'+h^2).^(3/2));

[U,S,V]=svd(Q);
s=diag(S);

M=0:maxM;
est=L * h.^(M+1)/(4*pi) .* ((2*M+1)/(1-h)+2/(1-h).^2)

reached=zeros(size(M));
K=zeros(size(X*Y'));
for M=0:maxM
  K=K+(2*M+1)/(4*pi)*h^M*my_legendre(M,X*Y');
  reached(M+1)=norm(Q-K);
end;

M=0:maxM;
semilogy(0:L-1,s/s(1),(M+1).^2,est/s(1),'x',(M+1).^2,reached/s(1),'+');

print temp.eps -deps
!ps2pdf temp.eps poisson_h_0.1.pdf 
!rm temp.eps

% bestapprox=zeros(L,1);
% for r=1:L
%   TS=diag([s(1:(r-1));zeros(L-r+1,1)]);
%   bestapprox(r)=norm(Q-U*TS*V');
% end;





