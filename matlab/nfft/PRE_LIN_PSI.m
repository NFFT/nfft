function f = PRE_LIN_PSI()
% If this flag is set, the convolution step (the multiplication with the sparse
% matrix B) uses linear interpolation from a lookup table of equispaced samples
% of the window function instead of exact values of the window function.
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

f = bitshift(1, 2);
