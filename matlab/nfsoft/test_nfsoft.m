N = 64;     % bandwidth (polynomial degree)
M = 100;	% number of points

% random nodes
alpha = rand(M,1)*2*pi;
beta = rand(M,1)*pi;
gamma = rand(M,1)*2*pi;
x = [alpha';beta';gamma'];

% random Fourier coefficients
fh = rand(nfsoft_f_hat_size(N),1);

plan = nfsoft_init(N,M);	% create plan
nfsoft_set_x(plan,x);       % set nodes
nfsoft_precompute(plan);	% node-dependent precomputation
nfsoft_set_f_hat(plan,fh);	% set Fourier coefficients
nfsoft_trafo(plan);         % transform
f = nfsoft_get_f(plan);     % get function values
nfsoft_finalize(plan);

%%
% nfsoftmex('set_f',plan,ones(M,1));
% nfsoftmex('precompute',plan);
% nfsoftmex('adjoint',plan);
% fh = nfsoftmex('get_fh',plan)
