function p = nfft_init_2d(N1,N2,M)
% Initialise plans
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

p = nfft('init_2d',N1,N2,M);
