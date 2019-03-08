% Example for the usage of class infft for the quadratic setting M=N.
%
% This is the testfile for the quadratic setting. For underdetermined and
% overdetermined cases see files 'test_underdetermined.m' and
% 'test_overdetermined.m' respectively.

clear

% Number of nodes
M = 1024;

% Random Fourier coefficients in [1,100]
fhat = rand(M,1)*99+1;

% Jittered nodes in [-0.5,0.5)
y = (-0.5:1/M:0.5-1/M) + 1/(4*M)*rand(1,M);

% Evaluations of a trigonometric polynomial at points y
f = exp(2*pi*1i*y'*(-M/2:M/2-1))*fhat;

%% Fast computation

% Initialization and node-dependent precomputations
% Using the default values. For customized settings see README.
plan = infft(y);
% Could also be computed as 
%   plan = infft(y,M);

plan.f = f; % Set function values
infft_trafo(plan); % Compute inverse nonequispaced Fourier transform

%% Direct computation

infft_trafo_direct(plan); % Compute samples directly

%% Output

% Graphical representation of the pointwise reconstruction
figure(1)
plot(-M/2:M/2-1,fhat,'o',-M/2:M/2-1,real(plan.fcheck),'*')
title('Pointwise reconstruction')
xlabel('$k$','Interpreter','latex')
ylabel('Fourier coefficients $\hat f_k$','Interpreter','latex')
legend('Fourier coefficients','approximation','Orientation','horizontal','Location','bestoutside')
xlim([-M/2-1,M/2])

% Graphical representation of pointwise errors
figure(2)
semilogy(-M/2:M/2-1,abs(plan.fcheck(:)-fhat(:)),'-sg',-M/2:M/2-1,abs(plan.fcheck(:)-fhat(:))./norm(fhat,Inf),'-dr')
title('Pointwise maximum errors')
xlabel('$k$','Interpreter','latex')
ylabel('pointwise errors')
legend('absolute maximum error','relative maximum error','Location','best')
xlim([-M/2-1,M/2])

% Computation of errors
err_abs_max = norm(plan.fcheck(:)-fhat(:),Inf);                  % Absolute ell_infinity error
err_rel_max = norm(plan.fcheck(:)-fhat(:),Inf)/norm(fhat,Inf);   % Relative ell_infinity error
err_abs_2 = norm(plan.fcheck(:)-fhat(:),2);                      % Absolute ell_2 error
err_rel_2 = norm(plan.fcheck(:)-fhat(:),2)/norm(fhat,2);         % Relative ell_2 error

% Output in command window
fprintf(['The absolute maximum error is ',num2str(err_abs_max,'%1.4e'),' and the relative maximum error is ',num2str(err_rel_max,'%1.4e'),'.\n\n'])