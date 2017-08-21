N = 32;     % bandwidth (polynomial degree)
M = nfsoft_f_hat_size(N);	% number of points

% random nodes
alpha = rand(M,1)*2*pi;
beta = acos(rand(M,1)*2-1);
gamma = rand(M,1)*2*pi;
x = [alpha';beta';gamma'];
M = length(alpha);

%% random Fourier coefficients
fh = rand(nfsoft_f_hat_size(N),1)./(1:nfsoft_f_hat_size(N)).';

plan = nfsoft_init(N,M,NFSOFT_NORMALIZED);  % create plan
nfsoft_set_x(plan,x);           % set nodes
nfsoft_precompute(plan);        % node-dependent precomputation
nfsoft_set_f_hat(plan,fh);      % set Fourier coefficients
nfsoft_trafo(plan);             % transform
f = nfsoft_get_f(plan);         % get function values
nfsoft_finalize(plan);

%% adjoint
plan = nfsoft_init(N,M,NFSOFT_NORMALIZED);
nfsoft_set_x(plan,x);
nfsoft_set_f(plan,f/M);
nfsoft_precompute(plan);
nfsoft_adjoint(plan);
fh2 = nfsoft_get_f_hat(plan);