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
N = 128;

% threshold (affects accuracy, 1000 is the default)
kappa = 1000;

% random Fourier coefficients
fh = f_hat(rand((N+1)*(N+1),1)./(1:(N+1)*(N+1))');

% precomputation
tic
nfsft_precompute(N,kappa);
fprintf('Time of precomputation:    %g seconds\n', toc);

% number of nodes
M = 5000;

% random nodes
ph=rand(1,M)*2*pi;
th=rand(1,M)*pi;
X=[ph;th];

fprintf('Performing a fast and direct NFSFT with bandwidth N = %d and M = %d nodes\n', N,M);

% Create plan.
plan = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT));

% Set nodes.
nfsft_set_x(plan,X);

% set Fourier coefficients
nfsft_set_f_hat(plan,double(fh));

% direct (slow) transform
tic
ndsft_trafo(plan);
fprintf('Time of direct transform:  %g seconds\n', toc);

% get function values
f_direct = nfsft_get_f(plan);

% fast NFSFT transform
tic
nfsft_trafo(plan);
fprintf('Time of fast transform:    %g seconds\n', toc);
% get function values
f_fast = nfsft_get_f(plan);


% random function values for adjoint
f = rand(M,1);

% adjoint transform
nfsft_set_f(plan,f)
tic
ndsft_adjoint(plan);
fprintf('Time of direct adjoint:    %g seconds\n', toc);

% get function values
fh2_direct = f_hat(nfsft_get_f_hat(plan));

% Fast adjoint transform
nfsft_set_f(plan,f)
tic
nfsft_adjoint(plan);
fprintf('Time of fast adjoint:      %g seconds\n', toc);
fh2_fast = f_hat(nfsft_get_f_hat(plan));

% finalize plan
nfsft_finalize(plan);


% Algorithm in Matlab
tic
f_mat = legendre(0,cos(X(2,:)),'norm').' / sqrt(2*pi) * fh(0,0);
fh2_mat = zeros(size(fh.f_hat));
fh2_mat(1)= legendre(0,cos(X(2,:)),'norm') / sqrt(2*pi) * f;
for n=1:N
    Lp = legendre(n,cos(X(2,:)),'norm') / sqrt(2*pi);
    f_mat=f_mat + (flipud(Lp).* exp((-n:0)'*1i*X(1,:))).'* fh(n,-n:0)...
     + (Lp(2:end,:).* exp((1:n)'*1i*X(1,:))).'* fh(n,1:n);
    
    fh2_mat(n^2+1:n^2+n+1)= (flipud(Lp).* exp(-(-n:0)'*1i*X(1,:)))* f;
    fh2_mat(n^2+n+2:n^2+2*n+1)= (Lp(2:end,:).* exp(-(1:n)'*1i*X(1,:)))* f;
end
fh2_mat = f_hat(fh2_mat);
fprintf('Time in Matlab:            %g seconds\n', toc);

% display error
fprintf('Relative error trafo fast vs matlab:     %g\n', norm(f_fast-f_mat)/norm(f_mat));
fprintf('Relative error trafo direct vs matlab:   %g\n', norm(f_direct-f_mat)/norm(f_mat));
fprintf('Relative error trafo fast vs direct:     %g\n', norm(f_fast-f_direct)/norm(f_direct));
fprintf('Relative error adjoint fast vs matlab:   %g\n', norm(fh2_fast-fh2_mat)/norm(fh2_mat));
fprintf('Relative error adjoint direct vs matlab: %g\n', norm(fh2_direct-fh2_mat)/norm(fh2_mat));
fprintf('Relative error adjoint fast vs direct:   %g\n', norm(fh2_fast-fh2_direct)/norm(fh2_direct));

% destroy all precomputations
% nfsft_forget()
