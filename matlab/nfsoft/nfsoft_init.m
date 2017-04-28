%Initialize NFSOFT plan
% N ... polynomial degree (bandwidth)
% M ... number of nodes
function plan = nfsoft_init(N, M, nfsoft_flags)
narginchk(2,3);
if( nargin<3)
    nfsoft_flags=0;
end
plan = nfsoftmex('init',N,M,nfsoft_flags);