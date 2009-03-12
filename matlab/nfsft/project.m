%
% Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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
% $Id$
function err = project()

kappa = 1000;
tmax = 0;
err = zeros(tmax+1,1);

nfsft_precompute(2^tmax,kappa);

ind = 1;
for N=2.^(0:tmax)
  [X,W] = gl(N);
  M = size(X,2);
  plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);
  nfsft_set_x(plan,X);
  nfsft_precompute_x(plan);
  f = ff(X).*W;
  nfsft_set_f(plan,f);
%  nfsft('display',plan);
  ndsft_adjoint(plan);
  fh = f_hat(nfsft_get_f_hat(plan));
  nfsft_finalize(plan);

  M = 20;
  X = repmat([2;1]*pi,1,M) .* rand(2,M);
  plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);
  nfsft_set_x(plan,X);
  nfsft_precompute_x(plan);
  nfsft_set_f_hat(plan,double(fh));
%  nfsft('display',plan);
  ndsft_trafo(plan);
  fa = nfsft_get_f(plan)';
  fe = ff(X);
  err(ind) = norm(fe-fa)/norm(fe);
  ind = ind + 1;
  nfsft_finalize(plan);
end

nfsft_forget();

% loglog(2.^(0:tmax),err);

% clear all;
