% Copyright (c) 2002, 2020 Jens Keiner, Stefan Kunis, Daniel Potts
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

% Test script of class nfct for spatial dimension d=1.
clear all;

M=19; % number of nodes
N=43; % number of Fourier coefficients in first direction

x=0.5*rand(M,1); %nodes

% Initialisation
plan=nfct(1,N,M); % create plan of class type nfct
%n=2^(ceil(log(N)/log(2))+1);
%plan=nfct(1,N,M,n,7,bitor(PRE_PHI_HUT,PRE_PSI),FFTW_MEASURE); % use of nfct_init_guru

plan.x=x; % set nodes in plan
%nfct_precompute_psi(plan); % precomputations

% NFCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhat=rand(N,1); % Fourier coefficients
fhatv=fhat(:);

% Compute samples with NFCT
plan.fhat=fhatv; % set Fourier coefficients
nfct_trafo(plan); % compute nonequispaced Fourier transform
f1=plan.f; % get samples

% Compute samples direct
k1=(0:N-1).';
f2=zeros(M,1);
for j=1:M
	x1j=x(j,1);
	f2(j)=sum( fhatv.*cos(2*pi*k1*x1j) );
end %for

% Compare results
max(abs(f1-f2))

% Adjoint NFCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation with NFCT
nfct_adjoint(plan);
fhat1=plan.fhat;

% Direct computation
fhat2=zeros(N,1);
for j=1:N
	k1j=k1(j);
	fhat2(j)=sum( plan.f.*cos(2*pi*k1j*x(:,1)) );
end %for

% Compare results
max(abs(fhat1-fhat2))

