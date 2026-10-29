classdef nfsoftTestcaseDelegateGLQuadratureOnline < nfsoftTestcaseDelegate
  %nfsoftTESTCASEDELEGATE Summary of this class goes here
  %   Detailed explanation goes here

% Copyright (c) 2018 The NFFT Authors
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <https://www.gnu.org/licenses/>.
  
  properties(Hidden=true,SetAccess='protected',GetAccess='public')
    f2 = [];
  end
  
  methods
    function h = nfsoftTestcaseDelegateGLQuadratureOnline(N)
      h.N = N;
    end
    
    function h = setup(h)
      fprintf('%-31s', 'nfsoft_gl_quadrature_online');
      N = h.N;
      
      [beta, W] = lgwt(N+1,-1,1);
      [alpha,beta,gamma]=meshgrid(0:2*N,beta,0:2*N);
      alpha=alpha(:)*2*pi/(2*N+1);
      beta =acos(beta(:));
      gamma=gamma(:)*2*pi/(2*N+1);
      X = [alpha';beta';gamma'];
      h.x = X;
      W = repmat(W,(2*N+1)^2,1)/2/(2*N+1)^2 * 8*pi^2;
      M = length(alpha);	% number of points
      h.M = M;
      
      fprintf(' N = %-5d,', h.N);
      fprintf(' M = %-5d,', h.M);
      
      fprintf(' nthreads = %d', nfsoft_get_num_threads);
      
      % random Fourier coefficients
      fh = rand(nfsoft_f_hat_size(N),1);
      h.f_hat = fh;
      
      % Create plan of class nfsoft.
      plan = nfsoft(N,M,bitor(NFSOFT_NORMALIZED,NFSOFT_USE_DPT));
      plan2 = nfsoft(N,M,NFSOFT_NORMALIZED);
      plan.x = X;
      plan2.x = X;

      plan.fhat = fh;
      plan2.fhat = fh;
          
      % nfsoft transform
      nfsoft_trafo(plan);
      nfsoft_trafo(plan2);
          
      % function values multiplied by quadrature weights
      h.f = plan.f .* W;
      h.f2 = plan2.f .* W;
    end
    
    function h = destroy(h)
      h.x = [];
      h.f_hat = [];
      h.f = [];
      h.f2 = [];
      h = destroy@nfsoftTestcaseDelegate(h);
    end
  end
  
end

