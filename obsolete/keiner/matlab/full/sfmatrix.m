function Y = sfmatrix(M,theta,phi,normalize)
% SFMATRIX Spherical Discrete Fourier transform matrix.
%   Y = SFMATRIX(M,THETA,PHI,NORMALIZE) computes the matrix corresponding to the 
%   evaluation of a bandlimited function with bandwidth M at the nodes 
%   (THETA,PHI). Y*A gives the function values at the nodes where A is a 
%   column vector containing the Fourier coefficients a_k^n, k = 0,...,M, 
%   n = -k,...,k in the order 
%
%     a_0^0,a_1^(-1),a_1^(0),a_1^(1),a_2^(-2),...,a_M^(M-1),a_M^(M).  
%
% The parameter NORMALIZE defines the normalization used for the spherical harmonics.
% Admissible values are 'unnorm', 'semi' and 'norm'.

  Y = zeros((M+1)*(M+1),length(theta));
  for k=0:M
    Y(k*k+1:k*k+2*k+1,:) = sheval(k,theta,phi,normalize);
  end
  Y = Y';
 
%   F = [];
%   [x,y] = meshgrid((-N):N,theta);
%   f1 = exp(-i*x.*y);
%   [x2,y2] = meshgrid((-N):N,phi);
%   f2 = exp(-i*x2.*y2);
%   for n=1:length(theta)
%     F = [F;kron(f2(n,:),f1(n,:))];   
%   end
  %exp(-2*pi*theta)
  