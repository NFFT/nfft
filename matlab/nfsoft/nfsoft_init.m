%Initialize NFSOFT plan
% plan = nfsoft_init(N, M)
% plan = nfsoft_init(N, M, nfsoft_flags)
% plan = nfsoft_init(N, M, nfsoft_flags, nfft_flags, nfft_cutoff, fpt_kappa)
% N ... polynomial degree (bandwidth)
% M ... number of nodes
% nfsoft_flags  ... can be NFSOFT_NORMALIZED, NFSOFT_REPRESENT,
% NFSOFT_USE_DPT, NFSOFT_USE_NDFT (default=0)

function plan = nfsoft_init(N, M, ...
    nfsoft_flags, nfft_flags, nfft_cutoff, fpt_kappa, fftw_size)
narginchk(2,7);
if(nargin<7 || isempty(fftw_size))
    fftw_size=4*N+4; % oversampling of NFFT
end
if(nargin<6 || isempty(fpt_kappa))
    fpt_kappa=1000;
end
if(nargin<5 || isempty(nfft_cutoff))
    nfft_cutoff=6;
end
if(nargin<4 || isempty(nfft_flags))
    nfft_flags=0;
end
if(nargin<3 || isempty(nfsoft_flags))
    nfsoft_flags=0;
end
plan = nfsoftmex('init',N,M,nfsoft_flags,nfft_flags,nfft_cutoff,...
    fpt_kappa, fftw_size);