% Example for the usage of class infft_adjoint.

clear

% Number of nodes
M = 1024;

% Number of Fourier coefficients
N = 2048;

% Random function values in [1,100]
f = rand(M,1)*99+1;

% Jittered nodes in [-0.5,0.5)
y = (-0.5:1/M:0.5-1/M) + 1/(4*M)*rand(1,M);

% Evaluations of a trigonometric polynomial at points y
fhat = f'*exp(-2*pi*1i*y'*(-N/2:N/2-1));

%% Fast computation

% Initialization and node-dependent precomputations
% Using the default values. For customized settings see README.
plan = infft(y,N);

plan.fhat = fhat; % Set Fourier coefficients
infft_adjoint(plan); % Compute inverse adjoint nonequispaced Fourier transform

%% Direct computation

infft_adjoint_direct(plan); % Compute samples directly

%% Output

% Graphical representation of the pointwise reconstruction
figure(1)
plot(y,f,'o',y,real(plan.ftilde),'*')
title('Pointwise reconstruction')
xlabel('$y$','Interpreter','latex')
ylabel('Evaluations at points $y$','Interpreter','latex')
legend('Function values','approximation','Orientation','horizontal','Location','bestoutside')
xlim([-0.5,0.5])

% Graphical representation of pointwise errors
figure(2)
semilogy(y,abs(plan.ftilde-f),'-sg',y,abs(plan.ftilde-f)./norm(f,Inf),'-dr')
title('Pointwise maximum errors')
xlabel('$y$','Interpreter','latex')
ylabel('pointwise errors')
legend('absolute maximum error','relative maximum error','Location','best')
xlim([-0.5,0.5])

% Computation of errors
err_abs_max = norm(plan.ftilde-f,Inf);               % Absolute ell_infinity error
err_rel_max = norm(plan.ftilde-f,Inf)/norm(f,Inf);   % Relative ell_infinity error
err_abs_2 = norm(plan.ftilde-f,2);                   % Absolute ell_2 error
err_rel_2 = norm(plan.ftilde-f,2)/norm(f,2);         % Relative ell_2 error

% Output in command window
fprintf(['The absolute maximum error is ',num2str(err_abs_max,'%1.4e'),' and the relative maximum error is ',num2str(err_rel_max,'%1.4e'),'.\n\n'])