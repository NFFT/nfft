
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

% Test script of class nfst for spatial dimension d=2.
clear all;

M=16; % number of nodes
N1=24; % number of Fourier coefficients in first direction
N2=32; % number of Fourier coefficients in second direction
N=[N1;N2];

x=0.5*rand(M,2); %nodes

% Initialisation
plan=nfst(2,N,M); % create plan of class type nfst
%n=2^(ceil(log(max(N))/log(2))+1);
%plan=nfst(2,N,M,n,n,7,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE); % use of nfst_init_guru

plan.x=x; % set nodes in plan
%nfst_precompute_psi(plan); % precomputations

% NFCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhat=rand(N1-1,N2-1); % Fourier coefficients
fhatv=fhat(:);

% Compute samples with NFCT
plan.fhat=fhatv; % set Fourier coefficients
nfst_trafo(plan); % compute nonequispaced Fourier transform
f1=plan.f; % get samples

% Compute samples direct
k1=1:N1-1;
k2=1:N2-1;
[K1,K2]=ndgrid(k1,k2);
k1=K1(:); clear K1;
k2=K2(:); clear K2;
f2=zeros(M,1);
for j=1:M
	x1j=x(j,1);
	x2j=x(j,2);
	f2(j)=sum( fhatv.*sin(2*pi*k1*x1j).*sin(2*pi*k2*x2j) );
end %for

% Compare results
max(abs(f1-f2))

% Adjoint NFCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation with NFCT
nfst_adjoint(plan);
fhat1=plan.fhat;

% Direct computation
fhat2=zeros((N1-1)*(N2-1),1);
for j=1:(N1-1)*(N2-1)
	k1j=k1(j);
	k2j=k2(j);
	fhat2(j)=sum( plan.f.*sin(2*pi*k1j*x(:,1)).*sin(2*pi*k2j*x(:,2)) );
end %for

% Compare results
max(abs(fhat1-fhat2))

