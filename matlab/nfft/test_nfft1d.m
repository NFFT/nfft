
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

% Test script of class nfft for spatial dimension d=1.
clear all;

M=16; % number of nodes
N=24; % number of Fourier coefficients in first direction

x=rand(M,1)-0.5; %nodes

% Initialisation
plan=nfft(1,N,M); % create plan of class type nfft
%n=2^(ceil(log(N)/log(2))+1);
%plan=nfft(1,N,M,n,8,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE); % use of nfft_init_guru

plan.x=x; % set nodes in plan and perform precomputations

% NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhat=rand(N,1); % Fourier coefficients
fhatv=fhat(:);

% Compute samples with NFFT
plan.fhat=fhatv; % set Fourier coefficients
nfft_trafo(plan); % compute nonequispaced Fourier transform
f1=plan.f; % get samples

% Compute samples direct
k1=(-N/2:N/2-1).';
f2=zeros(M,1);
for j=1:M
	x1j=x(j,1);
	f2(j)=sum( fhatv.*exp(-2*pi*1i*k1*x1j) );
end %for

% Compare results
max(abs(f1-f2))

% Adjoint NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation with NFFT
nfft_adjoint(plan);
fhat1=plan.fhat;

% Direct computation
fhat2=zeros(N,1);
for j=1:N
	k1j=k1(j);
	fhat2(j)=sum( plan.f.*exp(2*pi*1i*k1j*x(:,1)) );
end %for

% Compare results
max(abs(fhat1-fhat2))

