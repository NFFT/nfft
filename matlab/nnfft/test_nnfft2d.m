
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

% Test script of class nnfft for spatial dimension d=2.
clear all;

%M=16; % number of nodes
%N_1=24; % number of Fourier coefficients in first direction
%N_2=32; % number of Fourier coefficients in second direction
%N=[N_1;N_2];
%N_total=N_1*N_2; % total number of Fourier coefficients
M=4;
N_1=2;
N_2=2;
N=[N_1;N_2];
N_total=6;



x=rand(M,2)-0.5; %nodes
v=rand(N_total,2)-0.5; %nodes

% Plan initialisation simple interface
plan=nnfft(2,N_total,M,N); % create plan of class type nnfft


% Plan initialisation guru interface
%sigma=2; % oversampling factor
%N1_1=sigma*N_1; % FFTW length, must be even natural number!
%N1_2=sigma*N_2; % FFTW length, must be even natural number!
%m=6; % window cut-off parameter
%plan=nnfft(2,N_total,M,N,N1_1,N1_2,m,bitor(PRE_PHI_HUT,PRE_PSI)); % create plan of class type nnfft

plan.x=x; % set nodes in plan
plan.v=v; % set nodes in plan
nnfft_precompute_psi(plan); % precomputations

% NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fhat=rand(N_1,N_2); % Fourier coefficients
%fhat=rand(N_total,1);
%fhat=ones(N_total,1);

%Test mit zweitem Einheitsvektor
fhat=zeros(N_total,1);
fhat(2)=1;
fhatv=fhat;%(:);

% Compute samples with NNFFT
 plan.fhat=fhatv; % set Fourier coefficients
 nnfft_trafo(plan); % compute nonequispaced Fourier transform
 f1=plan.f; % get samples

% Compute samples direct
nnfft_trafo_direct(plan);
f2=plan.f; 

%% Compare results
disp('NNFFT vs NNDFT');
max(abs(f1-f2))


A=exp(-2*pi*1i*x*diag(N)*v');
f3=A*fhatv;

%tmpv=(v.*repmat(N',size(v,1),1)).';

%geflipter zweiter Einheitsvektor
%f3=exp(-2*pi*1i*(x*tmpv))*flipud(plan.fhat);
%f3=exp(-2*pi*1i*(tmpx*v.'))*fhatv;
max(abs(f2-f3))
