%TEST_NFSFT Example program:  of fast algorithm
%Use simple_test.m for the basic usage exampe.

% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%
fprintf('Number of threads: %d\n', nfsft_get_num_threads());

% maximum degree (bandwidth) of spherical harmonics expansions
N = 1024;
direct = 0;

% threshold (affects accuracy, 1000 is the default)
kappa = 10000;

% random Fourier coefficients
fh = f_hat(rand((N+1)*(N+1),1)./(1:(N+1)*(N+1))');

% precomputation
tic
nfsft_precompute(N,kappa,0,0*FPT_NO_FAST_ALGORITHM);
fprintf('Time of precomputation:    %g seconds\n', toc);

% number of nodes
M = (2*N+2)^2;

% nodes
ph=(-N-1:N)/(2*N+2)*2*pi;
th=(-N-1:N)/(2*N+2)*2*pi;
[ph,th]=meshgrid(ph,th);
X=[ph(:)';th(:)'];
index = th(:)>=0;

fprintf('Performing a fast and direct NFSFT with bandwidth N = %d and M = %d nodes\n', N,M);

% Create plan.
plan = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT));
plan2 = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT)+2^17);

% Set nodes.
nfsft_set_x(plan,X);
nfsft_set_x(plan2,X);

% set Fourier coefficients
nfsft_set_f_hat(plan,double(fh));
nfsft_set_f_hat(plan2,double(fh));

% direct (slow) transform
if(direct)
tic
ndsft_trafo(plan);
fprintf('Time of direct transform:  %g seconds\n', toc);
end

% get function values
f_direct = nfsft_get_f(plan);

% fast NFSFT transform
tic
nfsft_trafo(plan);
fprintf('Time of fast transform:    %g seconds\n', toc);
% get function values
f_fast = nfsft_get_f(plan);

% fast FSFT transform
tic
nfsft_trafo(plan2);
fprintf('Time of FSFT transform:    %g seconds\n', toc);
% get function values
f_fsft = nfsft_get_f(plan2);


% random function values for adjoint
f = rand(M,1);
f = f.*index;

% adjoint transform
if(direct)
nfsft_set_f(plan,f)
tic
ndsft_adjoint(plan);
fprintf('Time of direct adjoint:    %g seconds\n', toc);
end

% get function values
fh2_direct = f_hat(nfsft_get_f_hat(plan));

% Fast adjoint transform
nfsft_set_f(plan,f)
tic
nfsft_adjoint(plan);
fprintf('Time of fast adjoint:      %g seconds\n', toc);
fh2_fast = f_hat(nfsft_get_f_hat(plan));

% FSFT adjoint transform
nfsft_set_f(plan2,f)
tic
nfsft_adjoint(plan2);
fprintf('Time of FSFT adjoint:      %g seconds\n', toc);
fh2_fsft = f_hat(nfsft_get_f_hat(plan2));

% finalize plan
nfsft_finalize(plan);
nfsft_finalize(plan2);


% % Algorithm in Matlab
% tic
% f_mat = legendre(0,cos(X(2,:)),'norm').' / sqrt(2*pi) * fh(0,0);
% fh2_mat = zeros(size(fh.f_hat));
% %fh2_mat(1)= legendre(0,cos(X(2,:)),'norm') / sqrt(2*pi) * f;
% for n=1:N
%     Lp = legendre(n,cos(X(2,:)),'norm') / sqrt(2*pi);
%     f_mat=f_mat + (flipud(Lp).* exp((-n:0)'*1i*X(1,:))).'* fh(n,-n:0)...
%      + (Lp(2:end,:).* exp((1:n)'*1i*X(1,:))).'* fh(n,1:n);
%     
% %    fh2_mat(n^2+1:n^2+n+1)= (flipud(Lp).* exp(-(-n:0)'*1i*X(1,:)))* f;
% %    fh2_mat(n^2+n+2:n^2+2*n+1)= (Lp(2:end,:).* exp(-(1:n)'*1i*X(1,:)))* f;
% end
% fh2_mat = f_hat(fh2_mat);
% fprintf('Time in Matlab:            %g seconds\n', toc);

% display error
fprintf('Relative error trafo fast vs FSFT:     %g\n', norm(f_fast(index)-f_fsft(index))/norm(f_fast(index)));
% fprintf('Relative error trafo fast vs matlab:     %g\n', norm(f_fast(index)-f_mat(index))/norm(f_mat(index)));
% fprintf('Relative error trafo direct vs matlab:   %g\n', norm(f_direct(index)-f_mat(index))/norm(f_mat(index)));
if(direct)
fprintf('Relative error trafo FSFT vs direct:     %g\n', norm(f_fsft(index)-f_direct(index))/norm(f_direct(index)));
end
fprintf('Relative error adjoint fast vs FSFT:   %g\n', norm(fh2_fast-fh2_fsft)/norm(fh2_fast));
% fprintf('Relative error adjoint fast vs matlab:   %g\n', norm(fh2_fast-fh2_mat)/norm(fh2_mat));
% fprintf('Relative error adjoint direct vs matlab: %g\n', norm(fh2_direct-fh2_mat)/norm(fh2_mat));
if(direct)
fprintf('Relative error adjoint FSFT vs direct:   %g\n', norm(fh2_fsft-fh2_direct)/norm(fh2_direct));
end

% destroy all precomputations
% nfsft_forget()
