function [C,D] = DCT3matrix(n)
[X,Y] = meshgrid(0:n-1,0:n-1);
C = cos(X.*(2*Y+1)*pi/(2*n));
D = eye(n);
D(1,1) = 0.5;