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
if(nargin<7)
    fftw_size=4*N+4; % oversampling of NFFT
if(nargin<6)
    fpt_kappa=1000;
if(nargin<5)
    nfft_cutoff=6;
if(nargin<4)
    nfft_flags=0;
if( nargin<3)
    nfsoft_flags=0;
end
end
end
end
end
plan = nfsoftmex('init',N,M,nfsoft_flags,nfft_flags,nfft_cutoff,...
    fpt_kappa, fftw_size);