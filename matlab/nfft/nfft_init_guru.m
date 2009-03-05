function p = nfft_init_guru(varargin)
% Initialise plans
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

varargin{:}

  p = nfft('init_guru',varargin);
%p = nfft('init_guru',1,10,9,20,4,PRE_PHI_HUT,FFTW_MEASURE);
