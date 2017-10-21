%SIMPLE_TEST Example program: Basic usage principles
%
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts

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
disp(sprintf('Number of threads: %d\n', nnfft_get_num_threads()));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('A simple one dimensional example');

% maximum degree (bandwidth)
N = 8;
%N1=sigma*N;
N_total=3;
% number of nodes
M = 17;
 
% nodes
x=rand(1,M)-0.5;
%xi=N*(rand(1,N)-0.5);
v=rand(1,N_total)-0.5;

% Create plan.
%plan=nnfft_init(1,N_total,M,N);
plan=nnfft_init_guru(1,N_total,M,N,2*N,6,bitor(PRE_PSI,PRE_PHI_HUT));

% Set nodes.
nnfft_set_x(plan,x);
nnfft_set_v(plan,v);


% node-dependent precomputation
nnfft_precompute_psi(plan);

% Fourier coefficients
f_hat = rand(N_total,1)+i*rand(N_total,1);

% Set Fourier coefficients.
nnfft_set_f_hat(plan,double(f_hat));

% transform
nnfft_trafo(plan);

% function values
f = nnfft_get_f(plan)

% finalize plan
%nnfft_finalize(plan);


%%%%%%%%%%%%%%%%%%%%%

nnfft_trafo_direct(plan);
f2=nnfft_get_f(plan)
% finalize plan
nnfft_finalize(plan);
%%%%%%%%%%%%%%%%%%%%%%%

%A=exp(-2*pi*i*x'*(-N_total/2:N_total/2-1));
A=exp(-2*pi*1i*x'*N*v);
f3=A*f_hat

disp('NNFFT vs NNDFT');
error_vector=f-f2;
error_linfl1=norm(f-f2,inf)/norm(f_hat,1)

disp('NNFFT vs Direct');
error_vector=f-f3;
error_linfl1=norm(f-f3,inf)/norm(f_hat,1)


