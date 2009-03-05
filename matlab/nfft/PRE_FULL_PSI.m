function f = PRE_FULL_PSI()
% If this flag is set, the convolution step (the multiplication with the sparse
% matrix B) uses (2m+2)^dM precomputed values of the window function, in
% addition indices of source and target vectors are stored.
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

f = bitshift(1, 5);
