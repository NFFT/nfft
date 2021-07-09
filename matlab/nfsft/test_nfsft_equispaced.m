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
N = 32;
direct = 1;

% threshold (affects accuracy, 1000 is the default)
kappa = 10000;

% random Fourier coefficients
fh = f_hat(rand((N+1)*(N+1),1)./(1:(N+1)*(N+1))');

% precomputation
tic
nfsft_precompute(N,kappa,0,0*FPT_NO_FAST_ALGORITHM);
fprintf('Time of precomputation:    %g seconds\n', toc);

% number of nodes
M = (2*N+2) * (N+2);

% nodes
ph=(-N-1:N)/(2*N+2)*2*pi;
th=(0:N+1)/(2*N+2)*2*pi;
[ph,th]=meshgrid(ph,th);
X=[ph(:)';th(:)'];

fprintf('Performing a fast and direct NFSFT with bandwidth N = %d and M = %d nodes\n', N,M);

% Create plan.
plan = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT));
plan2 = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT)+NFSFT_EQUISPACED);

% Set nodes.
nfsft_set_x(plan,X);

% set Fourier coefficients
nfsft_set_f_hat(plan,double(fh));
nfsft_set_f_hat(plan2,double(fh));

% direct (slow) transform
if(direct)
tic
ndsft_trafo(plan);
fprintf('Time of direct transform:  %g seconds\n', toc);

% get function values
f_direct = nfsft_get_f(plan);
end

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

% adjoint transform
if(direct)
nfsft_set_f(plan,f)
tic
ndsft_adjoint(plan);
fprintf('Time of direct adjoint:    %g seconds\n', toc);

% get function values
fh2_direct = f_hat(nfsft_get_f_hat(plan));
end

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

% display error
fprintf('Relative error trafo fast vs FSFT:       %g\n', norm(f_fast-f_fsft)/norm(f_fast));
if(direct)
fprintf('Relative error trafo FSFT vs direct:     %g\n', norm(f_fsft-f_direct)/norm(f_direct));
end
fprintf('Relative error adjoint fast vs FSFT:     %g\n', norm(fh2_fast-fh2_fsft)/norm(fh2_fast));
if(direct)
fprintf('Relative error adjoint FSFT vs direct:   %g\n', norm(fh2_fsft-fh2_direct)/norm(fh2_direct));
end

% destroy all precomputations
% nfsft_forget()
