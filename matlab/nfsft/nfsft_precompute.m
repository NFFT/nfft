function nfsft_precompute(N,kappa,nfsft_flags,fpt_flags)
% Precompute
%
% Copyright 2007 Jens Keiner

if ~exist('nfsft_flags','var')
  nfsft_flags = 0;
end

if ~exist('fpt_flags','var')
  fpt_flags = 0;
end

nfsft('precompute',N,kappa,nfsft_flags,fpt_flags)
