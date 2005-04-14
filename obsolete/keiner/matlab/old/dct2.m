function y = dct2(x)

% dct2 - compute DCT-II transform
%
%   y = dct2(x);
%
%   Compute
%       y[j] = \sum{ x[k]*cos(pi/n*j*(k+1/2)) }
%
%   Copyright (c) 2003 Gabriel Peyré

n = length(x);
y(2:2:2*n,1) = x;
y = [y ; zeros(2*n,1)];
y = real( fft(y) );
y = y(1:n);