% Example for the usage of infft.

clear

% Number of knots
N = 1024;             

% Random integer fourier coefficients in [1,100]
fhat = ceil(rand(N,1)*100);                     

% Jittered knots in [-0.5,0.5)
y = (-0.5:1/N:0.5-1/N) + 1/(4*N)*rand(1,N); 

% Evaluations of a trigonometric polynomial
f = exp(2*pi*1i*y'*(-N/2:N/2-1))*fhat;          
 

%% Fast computation

plan = infft(y);
% Using the default values; could also be computed as
%   infft(y,'m',4,'p',4,'n',2*N,'eps_I',4*4/(2*N),'sigma',2);

plan.f = f;
infft_trafo(plan);

fbar = plan.fbar;
times = plan.times;

err_max_abs = max(abs(fhat-fbar));                  % maximum absolute error
err_max_rel = max(abs(fhat-fbar)./abs(fhat));       % maximum relative error
err_mean_abs = sum(abs(fhat-fbar))/N;               % mean absolute error
err_mean_rel = sum(abs(fhat-fbar)./abs(fhat))/N;	% mean relative error

%% Direct computation

infft_direct(plan,f);

fbar_direct = plan.fbar_direct;
times.t_direct = plan.times.t_direct;