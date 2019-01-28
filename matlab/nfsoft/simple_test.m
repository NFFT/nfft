%SIMPLE_TEST Example program: Basic usage principles
% Computes the sum
% f(alpha,beta,gamma) = sum_{l=0}^N sum_{m,n=-l}^l f_hat(l,m,n) D_l^{m,n}(alpha,beta,gamma)
% in terms of Wiegner-D functions D^{m,n}_l on a set of arbitary nodes
% (alpha,beta,gamma), which are Euler angles on the rotation group SO(3).

N = 32;     % bandwidth (polynomial degree)
M = nfsoft_f_hat_size(N);	% number of points

% random nodes (Euler angles)
alpha = rand(M,1)*2*pi;
beta = acos(rand(M,1)*2-1);
gamma = rand(M,1)*2*pi;
x = [alpha';beta';gamma'];
M = length(alpha);

%% random Fourier (Wigner) coefficients
fh = rand(nfsoft_f_hat_size(N),1)./(1:nfsoft_f_hat_size(N)).';

tic
plan = nfsoft(N,M,NFSOFT_NORMALIZED);  % create NFSOFT plan
plan.x = x;                     % set nodes
plan.fhat = fh;                 % set Fourier (Wigner) coefficients
nfsoft_trafo(plan);             % NFSOFT transform
f = plan.f;                     % get function values
toc