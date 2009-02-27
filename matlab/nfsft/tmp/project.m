function err = project()

kappa = 1000;
tmax = 8;
Nv = 2.^(0:tmax);
err = zeros(tmax+1,1);

nfsft_precompute(2^tmax,kappa);
 
ind = 1;
for N=Nv
  [X,W] = gl(N);
  M = size(X,2);
  plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);
  nfsft_set_x(plan,X);
  nfsft_precompute_x(plan);
  f = ff(X).*W;
  nfsft_set_f(plan,f);
%  nfsft('display',plan);
  nfsft_adjoint(plan);
  fh = f_hat(nfsft_get_f_hat(plan));
  figure
  imagesc(log10(abs(double(fh))))
  colorbar
  nfsft_finalize(plan);

  M = 20;
  X = repmat([2;1]*pi,1,M) .* rand(2,M);
  plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);
  nfsft_set_x(plan,X);
  nfsft_precompute_x(plan);
  nfsft_set_f_hat(plan,double(fh));
%  nfsft('display',plan);
  nfsft_trafo(plan);
  fa = nfsft_get_f(plan)';
  fe = ff(X);
  err(ind) = norm(fe-fa)/norm(fe);
  ind = ind + 1;
  nfsft_finalize(plan);
end

nfsft_forget();

figure
loglog(Nv,err);

end

function y = ff(x)

n = size(x,2);
y = ones(1,n);
ind = x(2,:) > pi/2;
y(ind) = 1./sqrt(1 + 3*cos(x(2,ind)).^2);

end

