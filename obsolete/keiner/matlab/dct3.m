function y = dct3(x)

% dct2 - compute DCT-III transform
%
%   y = dct3(x);
%
%   Compute
%       y[j] = 1/2*x[0] + \sum{ x[k]*cos(pi/n*k*(j+1/2)) }
%
%   Copyright (c) 2003 Gabriel Peyré

n = length(x);
y = [x ; zeros(3*n,1)];
y = real( fft(y) );
y = y(2:2:2*n) - x(1)/2;