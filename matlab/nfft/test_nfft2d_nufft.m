
% Copyright (c) 2002, 2019 Jens Keiner, Stefan Kunis, Daniel Potts
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

% Test script of class nfft for spatial dimension d=2.
% Comparison with Matlab nufft implementation (Matlab 2020a+)
clear all;

% Terminate the script to avoid errors when running in old Matlab or Octave.
% It would otherwise cause an error because nufft is not found in Matlab <=2019
% or in Octave, which would stop the test sctipt.
if ~exist('nufft')
  warning('Script stopped. The function "nufft" is only supperoted by Matlab 2020a or newer.')
  return
end

M=16000; % number of nodes
N1=128; % number of Fourier coefficients in first direction
N2=128; % number of Fourier coefficients in second direction
N=[N1;N2];

x=rand(M,2)-0.5; %nodes

% Initialisation
tic
plan=nfft(2,N,M); % create plan of class type nfft
%n=2^(ceil(log(max(N))/log(2))+1);
%plan=nfft(2,N,M,n,n,8,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE); % use of nfft_init_guru

plan.x=x; % set nodes in plan and perform precomputations
fprintf('Time NFFT pre    = %g sec\n',toc)

% NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhat=rand(N1,N2); % Fourier coefficients
fhatv=fhat(:);

% Compute samples with NFFT
plan.fhat=fhatv; % set Fourier coefficients
tic
nfft_trafo(plan); % compute nonequispaced Fourier transform
fprintf('Time NFFT trafo  = %g sec\n',toc)
f1=plan.f; % get samples

% Compute samples direct
k1=-N1/2:N1/2-1;
k2=-N2/2:N2/2-1;
[K1,K2]=ndgrid(k1,k2);
k=[K1(:) K2(:)];
clear K1 K2;
f2=zeros(M,1);
for j=1:M
	f2(j)=sum( fhatv.*exp(-2*pi*1i*(k*x(j,:).')) );
end %for
fprintf('Time direct      = %g sec\n',toc)

% Compare results
fprintf('Error nfft vs direct = %g\n',max(abs(f1-f2)))

tic
ndft_trafo(plan); % compute nonequispaced Fourier transform
fprintf('Time NFFT_direct = %g sec\n',toc)

%% Compare with Matlab nufft implementation (Matlab 2020a+)
[y1,y2] = meshgrid(-N2/2:(N2/2-1),-N1/2:(N1/2-1));
y = [y2(:) y1(:)];
tic
f3 = nufftn(fhat,y,x);
fprintf('Time nufft       = %g sec\n',toc)
fprintf('Error nufft = %g\n',max(abs(f1-f3)))

%% Adjoint NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation with NFFT
nfft_adjoint(plan);
fhat1=plan.fhat;

% Direct computation
fhat2=zeros(N1*N2,1);
for j=1:N1*N2
	fhat2(j)=sum( plan.f.*exp(2*pi*1i*(x*k(j,:).')) );
end %for

% Compare results
fprintf('Error nfft_adjoint = %g\n',max(abs(fhat1-fhat2)))
