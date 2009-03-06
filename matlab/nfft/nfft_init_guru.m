function p = nfft_init_guru(varargin)
% Initialise plans, no error handling
% Matlab might run into a segmentation violation for wrong parameters
%
% nfft_init_guru(d,N1,...,Nd,M,n1,...,nd,m,nfft_flags,fftw_flags)
%
% d            spatial dimension
% N1,...,Nd    bandwidths
% M            number of nodes
% n1,...,nd    fft lengths
% m            cut-off parameter
% nfft_flags   PRE_PHI_HUT | {FG_PSI, PRE_LIN_PSI, PRE_FG_PSI, PRE_PSI,
%	       PRE_FULL_PSI} | FFT_OUT_OF_PLACE
% fftw_flags   {FFTW_ESTIMATE, FFTW_MEASURE}
%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis

p = nfft('init_guru',varargin);
