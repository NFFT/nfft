%% Initialize parameters
d=2;
% rng(97)
kernel='multiquadric';
c = 1;%1/sqrt(N);
m = 4;
p = 3;
n = 156;
eps_I = .0625;
eps_B = .0625;

%% Initialize fastsum
% plan=fastsum_init_guru(d,N,M,n,m,p,kernel,c,eps_I,eps_B);
plan=fastsum_init(d,n,p,kernel,c,eps_I,eps_B);

%% set nodes
N=2024;
M=2024;
r = sqrt(rand(N,1))*(0.25-eps_B/2);
phi = rand(N,1)*2*pi;
x = [r.*cos(phi) r.*sin(phi)];
alpha = rand(N,1) + 1i*rand(N,1);
r = sqrt(rand(M,1))*(0.25-eps_B/2);
phi = rand(M,1)*2*pi;
y = [r.*cos(phi) r.*sin(phi)];

%% compute
fastsummex('set_source',plan,x,alpha,m)
fastsum_set_y(plan,y,7)

tic
fastsum_trafo_direct(plan)   % direct computation
f_dir=fastsum_get_f(plan);
t_direct=toc;
tic

tic
fastsum_trafo(plan)        % fast computation
t_trafo=toc;
f = fastsum_get_f(plan);

[f_old,f_dir_old]=fastsum(x,alpha,y,kernel,c,m,n,p,eps_I,eps_B);

% plot 
figure(1)
scatter(y(:,1),y(:,2),120,abs(f./f_dir-1),'.')
colorbar
figure(2)
scatter(y(:,1),y(:,2),120,abs(f_old./f_dir_old-1),'.')
colorbar

%% finalize 
fastsummex('finalize',plan)