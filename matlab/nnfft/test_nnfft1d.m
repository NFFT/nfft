
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
N_total=N; % total number of Fourier coefficients

x=rand(M,1)-0.5; %nodes
v=rand(N,1)-0.5; %nodes

% Plan initialisation simple interface
plan=nnfft(1,N_total,M,N); % create plan of class type nfft

% Plan initialisation guru interface
%sigma=2; % oversampling factor
%N1=sigma*N; % FFTW length, must be even natural number!
%m=6; % window cut-off parameter
%plan=nnfft(1,N,M,N,N1,m,bitor(PRE_PHI_HUT,PRE_PSI)); % create plan of class type nfft

plan.x=x; % set nodes in plan
plan.v=v;

nnfft_precompute_psi(plan); % precomputations

% NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhat=rand(N,1); % Fourier coefficients
fhatv=fhat(:);

% Compute samples with NFFT
plan.fhat=fhatv; % set Fourier coefficients
nnfft_trafo(plan); % compute nonequispaced Fourier transform (in space and time)
f1=plan.f; % get samples

% Compute samples direct
nnfft_trafo_direct(plan);
f2=plan.f; 

% Compare results
disp('NNFFT vs NNDFT');
max(abs(f1-f2))



