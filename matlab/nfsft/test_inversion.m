%TEST_INVERSION Example program: 
% We recover the spherical harmonics coefficients f_hat of a spherical
% function
% f(theta_d,phi_d) = sum_{k=0}^N sum_{n=-k}^k f_hat(k,n) Y_k^n(theta_d,phi_d)
% on a set of arbitrary nodes (theta_d,phi_d) for d=1..M, in spherical 
% coordinates 0 <= phi_d <2*pi and 0 <= theta_d <= pi.
% We use 2 different approaches: i) based on the inversion of the NFFT
% followed by a conversion to spherical harmonics coefficients, and ii) the
% inversion of the NFSFT.

% Copyright (c) 2002, 2021 Jens Keiner, Stefan Kunis, Daniel Potts
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

% maximum degree (bandwidth) of spherical harmonics expansions
N = 128;

% number of nodes
M = 4*N^2;

% random nodes on the sphere
phi = rand(1,M)*2*pi;
theta = rand(1,M)*pi;
X = [phi;theta];

% Create plan of class NFSFT.
plan = nfsft(N,M,NFSFT_NORMALIZED);

% Set nodes in spherical coordinates x = [phi; theta]
plan.x = X; 

% random sphrical harmonics coefficients
fh = f_hat((rand((N+1)*(N+1),1)-.5)./(((1:(N+1)^2)+1).').^.5*2);
% fh(1,0) = 1;
plan.fhat = fh;

% NFSFT transform
nfsft_trafo(plan);

% function values
f = plan.f;

%% inverse NFFT to compute Fourier coefficients on torus
addpath('../nfft')
iterations = 20;
% X2 = [X, [mod(X(1,:)+pi,2*pi); -X(2,:)]];
% f2 = [f;f];

% NFFT with same parameters as used internally in the NFSFT
NN = 2*[N N]+2;
plan_nfft = nfft(2,NN,M, 4*[N N], 6,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_ESTIMATE);
plan_nfft.x = X'/(2*pi); % nodes on [-.5 .5]

tic

% fh_nfft = lsqr(@(x, transp_flag) afun(transp_flag, x, p_nfft),f);

fh_nfft = zeros(NN); % initial guess
plan_nfft.fhat = fh_nfft(:);
plan_nfft.nfft_trafo;
rk = f - plan_nfft.f;
plan_nfft.f = rk;
plan_nfft.nfft_adjoint;
pk = plan_nfft.fhat;
sk = pk;
ak=[]; bk=[];
res_nfft = zeros(1,iterations);
for k=1:iterations    % Engl, Hanke, Neubauer: residual rk=dk
  plan_nfft.fhat = pk;
  plan_nfft.nfft_trafo;
  qk = plan_nfft.f;
  akm1 = ak;
  ak = (sk'*sk)/(qk'*qk);
  fh_nfft = (fh_nfft(:) + ak*pk);
%   fh_nfft = reshape(fh_nfft,NN);
%   fh_nfft = (fh_nfft + flip(fh_nfft,2).*(-1).^(-N-1:N).')/2;
%   fh_nfft = fh_nfft(:);
  rk = rk - ak*qk;
  plan_nfft.f = rk;
  plan_nfft.nfft_adjoint;
  nskm1 = (sk'*sk);
  sk = plan_nfft.fhat;
  bkm1=bk;
  bk = (sk'*sk)/nskm1;
  pk = sk + bk*pk;
  res_nfft(k) = sqrt((rk'*rk)/length(rk)); 
end

fprintf('iNFFT time:  %.3g seconds.\n',toc)

% Evaluate trigonometric polynomial at Clenshaw-Curtis nodes
[X_cc,W] = cc(N);
plan_nfft = nfft(2, 2*[N N]+2, length(X_cc), 4*[N N], 6,bitor(PRE_PHI_HUT,NFFT_OMP_BLOCKWISE_ADJOINT),FFTW_ESTIMATE);
plan_nfft.x = X_cc' / (2*pi); % nodes on [-.5 .5]
plan_nfft.fhat = fh_nfft;
plan_nfft.nfft_trafo;
f_nfft = plan_nfft.f;

% adjoint transform, using quadrature weights to recover spherical harmonics coefficients
plan_cc = nfsft(N, length(X_cc), NFSFT_NORMALIZED);
plan_cc.x = X_cc; 
plan_cc.f = f_nfft.*W';
nfsft_adjoint(plan_cc);

fh2 = plan_cc.fhat;
fprintf('Relative error of reconstructed f_hat: %.3g,   residual: %.3g\n',norm(fh2.f_hat-fh.f_hat)/norm(fh.f_hat),res_nfft(end));

%% CGNE for NFSFT
tic
fh_rec = zeros(size(fh.f_hat));
plan.fhat = fh_rec;
plan.nfsft_trafo;
rk = f - plan.f;
plan.f = rk;
plan.nfsft_adjoint;
pk = plan.fhat;
pk = pk.f_hat;
sk = pk;
ak=[]; bk=[];
res = zeros(1,iterations);
for k=1:iterations    % Engl, Hanke, Neubauer: residual rk=dk
  plan.fhat = pk;
  plan.nfsft_trafo;
  qk = plan.f;
  akm1 = ak;
  ak = (sk'*sk)/(qk'*qk);
  fh_rec = (fh_rec + ak*pk);
  rk = rk - ak*qk;
  plan.f = rk;
  plan.nfsft_adjoint;
  nskm1 = (sk'*sk);
  sk = plan.fhat;
  sk = sk.f_hat;
  bkm1=bk;
  bk = (sk'*sk)/nskm1;
  pk = sk + bk*pk;
  res(k) = sqrt((rk'*rk)/length(rk));
end
fprintf('iNFSFT time: %.3g seconds.\n',toc)
fprintf('Relative error of reconstructed f_hat: %.3g,   residual: %.3g\n',norm(fh_rec-fh.f_hat)/norm(fh.f_hat),res(end));

%% Plot
figure(1)
scatter(X(1,:),X(2,:),15e4/N^2,real(f),'filled')
colorbar
title('function values on random nodes on S²')

figure(2)
scatter(X_cc(1,:),X_cc(2,:),3e4/N^2,real(f_nfft),'filled')
colorbar
title('reconstructed function values on CC grid on T²')

figure(3)
fh_nfft = reshape(fh_nfft,NN);
imagesc(abs(fh_nfft));
colorbar
title('reconstructed Fourier coefficients on T²')

figure(4)
imagesc(abs(double(fh)));
colorbar
title('given (random) spherical harmonics coefficients')

% function y = afun(transp_flag, x, plan)
%   if strcmp(transp_flag, 'transp')
%     plan.f = x;
%     plan.nfft_adjoint;
%     y = plan.fhat;
%   elseif strcmp(transp_flag, 'notransp')
%     plan.fhat = x;
%     plan.nfft_trafo;
%     y = plan.f;
%   end
% end
