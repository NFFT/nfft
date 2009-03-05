function f = PRE_PSI()
% If this flag is set, the convolution step (the multiplication with the sparse
% matrix B uses (2m+2)dM precomputed values of the window function.
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

f = bitshift(1, 4);
