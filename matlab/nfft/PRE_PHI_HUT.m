function f = PRE_PHI_HUT()
% If this flag is set, the deconvolution step (the multiplication with the 
% diagonal matrix D) uses precomputed values of the Fourier transformed window
% function.
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

f = bitshift(1, 0);
