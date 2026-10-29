%TEST_NFSOFT Example program: Comparison with direct computation

% Copyright (c) 2018 The NFFT Authors
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <https://www.gnu.org/licenses/>.
N = 16;     % bandwidth (polynomial degree)
nfsoft_flags = 0;
nfft_flags = bitshift(1,12);      % NFFT_OMP_BLOCKWISE_ADJOINT (default)
nfft_cutoff = 6;
fpt_kappa = 1000;   % threshold (affects accuracy, 1000 is the default)

% random nodes (Euler angles)
M = 19;	% number of points
alpha = rand(M,1)*2*pi;
beta = acos(rand(M,1)*2-1);
gamma = rand(M,1)*2*pi;
x = [alpha';beta';gamma'];

% random Fourier (Wigner) coefficients
fh = rand(nfsoft_f_hat_size(N),1)./(1:nfsoft_f_hat_size(N)).';
% fh = [1;1;0;0;0;0;0;0;0;0];

fftw_size = 4*N+4;
fprintf('Perform NFSOFT with bandwidth N = %d, M = %d nodes and oversampling %d\n',N,M,fftw_size)

%% create NFSOFT plan
tic
plan = nfsoft_init(N,M,nfsoft_flags,nfft_flags,nfft_cutoff,fpt_kappa,fftw_size);
nfsoft_set_x(plan,x);           % set nodes
fprintf('Time of initialization: %g seconds\n', toc);

nfsoft_set_f_hat(plan,fh);      % set Fourier (Wigner) coefficients

%% NFSOFT transform
tic
nfsoft_trafo(plan);
fprintf('Time of transform:      %g seconds\n', toc);

%% get function values
f = nfsoft_get_f(plan);

% adjoint transform
tic
nfsoft_adjoint(plan);
fprintf('Time of adjoint:        %g seconds\n', toc);
fh2 = nfsoft_get_f_hat(plan);
nfsoft_finalize(plan)

%% Direct transform in Matlab
tic
f_mat = zeros(size(f));
fh2_mat = zeros(size(fh));
for m=1:M
    for n=0:N
        D = wignerD(n,alpha(m),beta(m),gamma(m));
        D = D(:);
        f_mat(m) = f_mat(m) + D.' * fh(nfsoft_f_hat_size(n-1)+1:nfsoft_f_hat_size(n));
        fh2_mat(nfsoft_f_hat_size(n-1)+1:nfsoft_f_hat_size(n)) = ...
          fh2_mat(nfsoft_f_hat_size(n-1)+1:nfsoft_f_hat_size(n)) + conj(D) * f(m);
    end
end
fprintf('Time of matlab trafo:   %g seconds\n', toc);
fprintf('Relative error of trafo:   %g\n',norm(f_mat-f)/norm(f_mat));
fprintf('Relative error of adjoint: %g\n',norm(fh2_mat(:)-fh2(:))/norm(fh2_mat(:)));
