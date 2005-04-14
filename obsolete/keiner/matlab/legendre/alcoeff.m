function [alpha,beta,gamma] = alcoeff(k,n)
% ALCOEFF Three-term recurrence coefficients for associated Legendre functions
%   [ALPHA,BETA,GAMMA] = ALCOEFF(K,N) computes the three-term recurrence
%   coefficients \alpha_k^n, \beta_k^n and \gamma_k^n (k,n >= 0).

  % Parameter checking
  if (exist('k','var') == 0)
    error('alcoeff: Input parameter ''k'' undefined.');
  end
  if (exist('n','var') == 0)
    error('alcoeff: Input parameter ''n'' undefined.');
  end
  if (k < 0 || n < 0)
    error('alcoeff: Input parameters ''k'' and ''n'' must be greater or equal than zero.');
  end
  
  % Compute alpha
  if k < n
    alpha = (-1)^(k+1);
  else
    alpha = (2*k+1)/(sqrt((k-n+1)*(k+n+1)));
  end
  
  % Compute beta
  if k < n
    beta = 1;
  else
    beta = 0;
  end
  
  % Compute gamma
  if k <= n
    gamma = 0;
  else
    gamma = -(sqrt((k-n)*(k+n)))/(sqrt((k-n+1)*(k+n+1)));
  end