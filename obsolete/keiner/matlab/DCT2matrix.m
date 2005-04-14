function [C] = DCT2matrix(n)
[X,Y] = meshgrid(0:n-1,0:n-1);
C = cos(Y.*(2*X+1)*pi/(2*n));
