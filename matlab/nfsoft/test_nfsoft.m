%TEST_NFSOFT Example program: Advanced usage principles
N = 32;     % bandwidth (polynomial degree)
nfsoft_flags = NFSOFT_NORMALIZED; % L2-normalized Wigner functions (rotational harmonics)
% nfsoft_flags = bitor(NFSOFT_NORMALIZED,NFSOFT_USE_DPT);
nfft_flags = bitshift(1,12);      % NFFT_OMP_BLOCKWISE_ADJOINT (default)
nfft_cutoff = 6;
fpt_kappa = 1000;   % threshold (affects accuracy, 1000 is the default)

%% Gauss-Legendre quadrature nodes on SO(3)
% These allow a rectonstruction with the adjoint transform.
oldpath = addpath('../nfsft');
[beta, w] = lgwt(N+1,-1,1);
[alpha,beta,gamma]=meshgrid(0:2*N,beta,0:2*N);
alpha=alpha(:)*2*pi/(2*N+1);
beta =acos(beta(:));
gamma=gamma(:)*2*pi/(2*N+1);
x = [alpha';beta';gamma'];
w = repmat(w,(2*N+1)^2,1)/2/(2*N+1)^2 * 8*pi^2;
M = length(alpha);	% number of points

%% random Fourier (Wigner) coefficients
fh = rand(nfsoft_f_hat_size(N),1)./(1:nfsoft_f_hat_size(N)).';

% compare different oversamplings fftw_size (must be >= 2*N+2 and even, default 4*N+4)
% (higher oversampling means better accuracy but more time and memory required)
for fftw_size = [3*N, 4*N+4]
    fprintf('\nPerform NFSOFT with bandwidth N = %d and oversampling %d\n',N,fftw_size)

    % create NFSOFT plan
    tic
    plan = nfsoft(N,M,nfsoft_flags,nfft_flags,nfft_cutoff,fpt_kappa,fftw_size);
    plan.x = x;           % set nodes
    fprintf('Time of initialization: %g seconds\n', toc);

    plan.fhat = fh;      % set Fourier (Wigner) coefficients

    % NFSOFT transform
    tic
    nfsoft_trafo(plan);
    fprintf('Time of transform:      %g seconds\n', toc);

    % get function values
    f = plan.f;

    % adjoint transform with quadrature
    plan.f = f.*w;
    tic
    nfsoft_adjoint(plan);
    fprintf('Time of adjoint:        %g seconds\n', toc);

    fprintf('Relative error of reconstructed coefficients: %g\n',norm(plan.fhat-fh)/norm(fh));
end

path(oldpath);
