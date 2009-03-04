function p = nfft_init_1d(N,M)
% Initialise plans
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

p = nfft('init_1d',N,M);
