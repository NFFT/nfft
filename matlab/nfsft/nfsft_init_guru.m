function p = nfsft_init_guru(N,M,flags,nfft_flags,nfft_cutoff)
% Initialise plans (guru)
%
% Copyright (c) 2007 Jens Keiner

p = nfsft('init_guru',N,M,flags,nfft_flags,nfft_cutoff);
