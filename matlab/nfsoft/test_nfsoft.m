%TEST_NFSOFT Example program: Advanced usage principles
N = 32;     % bandwidth (polynomial degree)
nfsoft_flags = NFSOFT_NORMALIZED; % L2-normalized Wigner functions (rotational harmonics)
nfft_flags = 0;
nfft_cutoff = 6;
fpt_kappa = 1000;   % threshold (affects accuracy, 1000 is the default)

%% Gauss-Legendre quadrature nodes
% These allow a rectonstruction with the adjoint transform.
oldpath = addpath('../nfsft');
[beta, wb] = lgwt(N+1,-1,1);
[alpha,beta,gamma]=meshgrid(0:2*N,beta,0:2*N);
alpha=alpha(:)*2*pi/(2*N+1);
beta =acos(beta(:));
gamma=gamma(:)*2*pi/(2*N+1);
x = [alpha';beta';gamma'];
w = repmat(wb,(2*N+1)^2,1)/2/(2*N+1)^2 * 8*pi^2;
M = length(alpha);	% number of points

%% random Fourier (Wigner) coefficients
fh = rand(nfsoft_f_hat_size(N),1)./(1:nfsoft_f_hat_size(N)).';

% compare different oversamplings (must be >= 2*N+2 and even, default 8*N)
% (higher oversampling means better accuracy but more time and memory required)
for N_oversampled = [8*N, 4*N]
    fprintf('\nPerform NFSOFT with bandwidth N = %d and oversampling %d\n',N,N_oversampled)

    % create NFSOFT plan
    tic
    plan = nfsoft_init(N,M,nfsoft_flags,nfft_flags,nfft_cutoff,fpt_kappa,N_oversampled);
    nfsoft_set_x(plan,x);           % set nodes
    nfsoft_precompute(plan);        % node-dependent precomputation
    fprintf('Time of initialization: %g seconds\n', toc);

    nfsoft_set_f_hat(plan,fh);      % set Fourier (Wigner) coefficients

    % NFSOFT transform
    tic
    nfsoft_trafo(plan);
    fprintf('Time of transform:      %g seconds\n', toc);

    % get function values
    f = nfsoft_get_f(plan);

    % adjoint transform with quadrature
    nfsoft_set_f(plan,f.*w)
    tic
    nfsoft_adjoint(plan);
    fprintf('Time of adjoint:        %g seconds\n', toc);
    fh2 = nfsoft_get_f_hat(plan);

    % finalize plan
    nfsoft_finalize(plan)

    fprintf('Relative error of reconstructed coefficients: %g\n',norm(fh2-fh)/norm(fh));
end

path(oldpath);
