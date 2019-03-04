% Example for the usage of class infft for the underdetermined setting M<N.

clear

% Number of nodes
M = 1024;

% Number of Fourier coefficients
N = 2048;

% Jittered nodes in [-0.5,0.5)
y = (-0.5:1/M:0.5-1/M) + 1/(4*M)*rand(1,M);

% Evaluate trigonometric function at nonequispaced nodes
f = ((cos(pi*y.^2)).^2.*sin(10*y.^2))';

% Exactly computed Fourier coefficients
load('coeffs_2048.mat')

%% Fast computation

% Initialization and node-dependent precomputations
% Using the default values. For customized settings see README.
plan = infft(y,N);

plan.f = f; % Set function values
infft_trafo(plan); % Compute inverse nonequispaced Fourier transform

%% Direct computation

infft_direct(plan); % Compute samples directly

%% Output

% Graphical representation of the pointwise reconstruction
figure(1)
plot(-N/2:N/2-1,fhat,'o',-N/2:N/2-1,real(plan.fcheck),'*')
title('Pointwise reconstruction')
xlabel('$k$','Interpreter','latex')
ylabel('Fourier coefficients $\hat f_k$','Interpreter','latex')
legend('Fourier coefficients','approximation','Orientation','horizontal','Location','bestoutside')
xlim([-N/2-1,N/2])

% Graphical representation of pointwise errors
figure(2)
semilogy(-N/2:N/2-1,abs(plan.fcheck(:)-fhat(:)),'-sg',-N/2:N/2-1,abs(plan.fcheck(:)-fhat(:))./norm(fhat,Inf),'-dr')
title('Pointwise maximum errors')
xlabel('$k$','Interpreter','latex')
ylabel('pointwise errors')
legend('absolute maximum error','relative maximum error','Location','best')
xlim([-N/2-1,N/2])

% Computation of errors
err_abs_max = norm(plan.fcheck(:)-fhat(:),Inf);                 % Absolute ell_infinity error
err_rel_max = norm(plan.fcheck(:)-fhat(:),Inf)/norm(fhat,Inf);	% Relative ell_infinity error
err_abs_2 = norm(plan.fcheck(:)-fhat(:),2);                     % Absolute ell_2 error
err_rel_2 = norm(plan.fcheck(:)-fhat(:),2)/norm(fhat,2);        % Relative ell_2 error

% Output in command window
fprintf(['The absolute maximum error is ',num2str(err_abs_max,'%1.4e'),' and the relative maximum error is ',num2str(err_rel_max,'%1.4e'),'.\n\n'])