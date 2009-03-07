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