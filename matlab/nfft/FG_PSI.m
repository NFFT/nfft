function f = FG_PSI()
% If this flag is set, the convolution step (the multiplication with the
% sparse matrix B) uses particular properties of the Gaussian window function
% to trade multiplications for direct calls to exponential function.
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

f = bitshift(1, 1);
