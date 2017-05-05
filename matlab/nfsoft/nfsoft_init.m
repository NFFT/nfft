%Initialize NFSOFT plan
% plan = nfsoft_init(N, M)
% plan = nfsoft_init(N, M, nfsoft_flags)
% plan = nfsoft_init(N, M, nfsoft_flags, nfft_flags, nfft_cutoff, fpt_kappa)
% N ... polynomial degree (bandwidth)
% M ... number of nodes

function plan = nfsoft_init(N, M, nfsoft_flags, nfft_flags, nfft_cutoff, fpt_kappa)
narginchk(2,6);
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
plan = nfsoftmex('init',N,M,nfsoft_flags,nfft_flags,nfft_cutoff,fpt_kappa);