function err = projection()
% threshold
kappa = 1000;
% polynomial degrees
emin = 0; emax = 9; Nv = 2.^(emin:emax);
% saves errors
err = [Nv',zeros(length(Nv),1)];
% uniform random nodes
Me = 1000; xe = [2*pi*rand(1,Me);acos(2*rand(1,Me)-1)]; fe = ff(xe);
% precomputation
nfsft_precompute(max(Nv),kappa);

% loop over polynomial degrees
j = 1;
for N = Nv
  % projection using Gauss-Legendre points
  [x,w] = gl(N);
  M = size(x,2);
  plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);
  nfsft_set_x(plan,x);
  nfsft_precompute_x(plan);
  f = ff(x).*w;
  nfsft_set_f(plan,f);
  nfsft_adjoint(plan);
  fh = f_hat(nfsft_get_f_hat(plan));
  nfsft_finalize(plan);

  % evaluation at random nodes
  plan = nfsft_init_advanced(N,Me,NFSFT_NORMALIZED);
  nfsft_set_x(plan,xe);
  nfsft_precompute_x(plan);
  nfsft_set_f_hat(plan,double(fh));
  nfsft_trafo(plan);
  fa = nfsft_get_f(plan)';
  err(j,2) = norm(fe-fa)/norm(fe);
  j = j + 1;
  nfsft_finalize(plan);
end

% delete precomputed data
nfsft_forget();

% error plot
figure;
loglog(Nv,err(:,2));

end

% the function f
function y = ff(x)
n = size(x,2);
y = ones(1,n);
j = x(2,:) > pi/2;
y(j) = 1./sqrt(1 + 3*cos(x(2,j)).^2);
end

