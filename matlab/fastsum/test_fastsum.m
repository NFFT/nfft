%% Initialize parameters
d=2;
N=2000;
M=2000;
kernel='gaussian';
c = 1/sqrt(N);
m = 4;
p = 3;
n = 156;
eps_I = .0625;
eps_B = .0625;
r = sqrt(rand(N,1))*(0.25-eps_B/2);
phi = rand(N,1)*2*pi;
x = [r.*cos(phi) r.*sin(phi)];
alpha = rand(N,1) + 1i*rand(N,1);
r = sqrt(rand(M,1))*(0.25-eps_B/2);
phi = rand(M,1)*2*pi;
y = [r.*cos(phi) r.*sin(phi)];

%% Perform fastsum
plan=fastsum_init_guru(d,N,M,n,m,p,kernel,c,eps_I,eps_B);
fastsum_set_x(plan,x)
fastsum_set_y(plan,y)
fastsum_set_alpha(plan,alpha)

% fastsum_trafo_direct(plan)   % direct computation
% f_dir=fastsum_get_f(plan);

fastsum_precompute(plan)
fastsum_trafo(plan)        % fast computation
f = fastsum_get_f(plan);
b = fastsum_get_b(plan);   % get expansion coefficients
fastsum_finalize(plan)

%% plot 
figure(1)
scatter(x(:,1),x(:,2),120,real(alpha),'.')
colorbar
figure(2)
scatter(y(:,1),y(:,2),120,real(f),'.')
colorbar