function f_hat = nfft_get_f_hat(p)
% Get Fourier coefficients from plan
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

f_hat = nfft('get_f_hat',p);
