function y = aleval(k,n,c,x)
% ALEVAL Evaluate associated Legendre polynomials
%   Y = ALEVAL(K,N,C,X) evaluates the associated Legendre polynomial
%   P_K^N(\cdot,C) at the nodes X.

  % Parameter checking
  if (exist('k','var') == 0)
    error('aleval: Input parameter ''k'' undefined.');
  end
  if (exist('n','var') == 0)
    error('aleval: Input parameter ''n'' undefined.');
  end
  if (exist('c','var') == 0)
    error('aleval: Input parameter ''c'' undefined.');
  end
  if (exist('x','var') == 0)
    error('aleval: Input parameter ''x'' undefined.');
  end
  if (k < 0 || n < 0 || c < 0)
    error('aleval: Input parameters ''k'', ''n'' and ''c'' must be greater or equal than zero.');
  end

  % Compute the function value.
  if k == -1
    % P_{-1}^n(\cdot,C) = 0 
    y = zeros(size(x));
  else
    if k == 0
      % P_{-1}^n(\cdot,C) = 1 
      y = ones(size(x));
    else
      % P_{k+1}^n(\cdot,C) = (alpha_{k+c}^n * x + beta_{k+c}^n) *
      % P_{k}^n(\cdot,C) + gamma_{k+c}^n * P_{k-1}^n(\cdot,C)
      % Clenshaw
      a = ones(size(x));
      b = zeros(size(x));
      for j = k:-1:2
        a_old = a;
        [alpha,beta,gamma] = alcoeff(j-1+c,n);
        a = b + alpha*x.*a_old + beta*a_old;
        b = gamma*a_old;
      end    
      [alpha,beta,gamma] = alcoeff(c,n);
      y = alpha*x.*a + beta*a + b;
    end
  end