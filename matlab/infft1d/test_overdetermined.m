% Example for the usage of class infft for the overdetermined setting M>N.

clear

% Number of nodes
M = 2048;

% Number of Fourier coefficients
N = 1024;

% Random Fourier coefficients in [1,100]
fhat = rand(N,1)*99+1;

% Jittered nodes in [-0.5,0.5)
y = (-0.5:1/M:0.5-1/M) + 1/(4*M)*rand(1,M);

% Evaluations of a trigonometric polynomial at points y
f = exp(2*pi*1i*y'*(-N/2:N/2-1))*fhat;

%% Fast computation

% Initialization and node-dependent precomputations
% Using the default values. For customized settings see README.
plan = infft(y,N);

plan.f = f; % Set function values
infft_trafo(plan); % Compute inverse nonequispaced Fourier transform

%% Direct computation

infft_direct(plan); % Compute samples directly

%% Fast computation via Gohberg-Semencul formula (only for M>=N)

plan2 = infft(y,N,'flag_toeplitz',1);

plan2.f = f; % Set function values
infft_trafo(plan2); % Compute inverse nonequispaced Fourier transform

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
% Matrix approach
err_abs_max = norm(plan.fcheck(:)-fhat(:),Inf);                        % Absolute ell_infinity error
err_rel_max = norm(plan.fcheck(:)-fhat(:),Inf)/norm(fhat,Inf);         % Relative ell_infinity error
err_abs_2 = norm(plan.fcheck(:)-fhat(:),2);                            % Absolute ell_2 error
err_rel_2 = norm(plan.fcheck(:)-fhat(:),2)/norm(fhat,2);               % Relative ell_2 error
% Toeplitz approach
err_abs_max_toep = norm(plan2.fcheck(:)-fhat(:),Inf);                  % Absolute ell_infinity error
err_rel_max_toep = norm(plan2.fcheck(:)-fhat(:),Inf)/norm(fhat,Inf);   % Relative ell_infinity error
err_abs_2_toep = norm(plan2.fcheck(:)-fhat(:),2);                      % Absolute ell_2 error
err_rel_2_toep = norm(plan2.fcheck(:)-fhat(:),2)/norm(fhat,2);         % Relative ell_2 error

% Output in command window
fprintf(['Using the matrix approach the absolute maximum error is ',num2str(err_abs_max,'%1.4e'),' and the relative maximum error is ',num2str(err_rel_max,'%1.4e'),'.\n'])
fprintf(['Using the Gohberg-Semencul formula the absolute maximum error is ',num2str(err_abs_max_toep,'%1.4e'),' and the relative maximum error is ',num2str(err_rel_max_toep,'%1.4e'),'.\n\n'])