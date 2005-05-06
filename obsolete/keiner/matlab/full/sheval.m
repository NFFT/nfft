function y = sheval(k,theta,phi,normalize)
% SHEVAL Evaluate spherical harmonics
%   Y = SHEVAL(K,THETA,PHI,NORMALIZE) evaluates the spherical harmonics Y_K^N
%   for N = -K,..,K at the nodes (THETA,PHI), where 0 <= THETA <= pi, 
%   0 <= PHI <= 2*pi. 
%   The i-th row contains the function values of Y_K^{i-K-1}.
%   In general, the spherical harmonics Y_K^N are defined by
%
%     Y_K^N(THETA,PHI) = P_K^{|N|}(cos(THETA)) * e^{i*N*PHI}
%
%   where
%
%     P_K^N = \delta * (1-x^2)^(N/2) * (d/dX)^N * P_K,
%     P_K = (1/(2^K * K!)) * (d/dX)^K * (x^2 - 1)^K.
%
%   The parameter NORMALIZE sets the normalization used:
%
%     \delta = (-1)^N                          if NORMALIZE == 'unnorm',
%     \delta = ((K-N)!/(K+N)!)^(1/2)           if NORMALIZE == 'semi',
%     \delta = ((K+1/2) * (K-N)!/(K+N)!)^(1/2) if NORMALIZE == 'norm'.

  % Parameter checking
  if (exist('k','var') == 0)
    error('sheval: Input parameter ''k'' undefined.');
  end
  if (exist('theta','var') == 0)
    error('sheval: Input parameter ''theta'' undefined.');
  end
  if (exist('phi','var') == 0)
    error('sheval: Input parameter ''phi'' undefined.');
  end
  if (exist('normalize','var') == 0)
    error('sheval: Input parameter ''normalize'' undefined.');
  end
  
  % Compute matrix containing row-wise the function values P_K^N for N =
  % 0,...,K at the nodes THETA.
  if strcmp(normalize,'unnorm') == 1
    p = legendre(k,cos(theta),'unnorm');    
    enorm = 1.0;
  else
    if strcmp(normalize,'semi') == 1
      p = legendre(k,cos(theta),'sch');     
      if k>0
        p(2:end,:) = (1/sqrt(2))*p(2:end,:);
      end
      enorm = 1.0;
    else
      if strcmp(normalize,'norm') == 1
        p = legendre(k,cos(theta),'norm');     
        enorm = 1.0;
    else
        error('sheval: Wrong value for parameter ''normalize''.');
      end
    end
  end

  % Compute matrix containing row-wise the function values Y_K^N for
  % N = -K,...,K at the nodes (THETA,PHI).
  [A,B] = meshgrid(phi,(-k:1:k));
  y = [flipud(p(2:end,:));p] .* (enorm*exp(i*A.*B));
 
