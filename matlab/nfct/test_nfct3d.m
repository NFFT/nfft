
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

% Test script of class nfct for spatial dimension d=3.
clear all;

M=16; % number of nodes
N1=24; % number of Fourier coefficients in first direction
N2=32; % number of Fourier coefficients in second direction
N3=30; % number of Fourier coefficients in third direction
N=[N1;N2;N3];

x=0.5*rand(M,3); %nodes

% Initialisation
plan=nfct(3,N,M); % create plan of class type nfct
%n=2^(ceil(log(max(N))/log(2))+1);
%plan=nfct(3,N,M,n,n,n,7,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE); % use of nfct_init_guru

plan.x=x; % set nodes in plan
%nfct_precompute_psi(plan); % precomputations

% NFCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhat=rand(N1,N2,N3); % Fourier coefficients
fhatv=fhat(:);

% Compute samples with NFCT
plan.fhat=fhatv; % set Fourier coefficients
nfct_trafo(plan); % compute nonequispaced Fourier transform
f1=plan.f; % get samples

% Compute samples direct
k1=0:N1-1;
k2=0:N2-1;
k3=0:N3-1;
[K1,K2,K3]=ndgrid(k1,k2,k3);
k1=K1(:); clear K1;
k2=K2(:); clear K2;
k3=K3(:); clear K3;
f2=zeros(M,1);
for j=1:M
	x1j=x(j,1);
	x2j=x(j,2);
	x3j=x(j,3);
	f2(j)=sum( fhatv.*cos(2*pi*k1*x1j).*cos(2*pi*k2*x2j).*cos(2*pi*k3*x3j) );
end %for

% Compare results
max(abs(f1-f2))

% Adjoint NFCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation with NFCT
nfct_adjoint(plan);
fhat1=plan.fhat;

% Direct computation
fhat2=zeros(N1*N2*N3,1);
for j=1:N1*N2*N3
	k1j=k1(j);
	k2j=k2(j);
	k3j=k3(j);
	fhat2(j)=sum( plan.f.*cos(2*pi*k1j*x(:,1)).*cos(2*pi*k2j*x(:,2)).*cos(2*pi*k3j*x(:,3)) );
end %for

% Compare results
max(abs(fhat1-fhat2))

