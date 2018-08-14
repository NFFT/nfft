classdef nfsftTestcaseDelegateGLQuadratureOnline < nfsftTestcaseDelegate
  %NFSFTTESTCASEDELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Hidden=true,SetAccess='protected',GetAccess='public')
    f2 = [];
  end
  
  methods
    function h = nfsftTestcaseDelegateGLQuadratureOnline(N)
      h.N = N;
    end
    
    function h = setup(h)
      fprintf('%-31s', 'nfsft_gl_quadrature_online');
      N = h.N;

      [X,W] = gl(N);
      h.x = X;
      M = size(X,2);
      h.M = M;

      fprintf(' N = %-5d,', h.N);
      fprintf(' M = %-5d,', h.M);
      
      fprintf(' nthreads = %d', nfsft_get_num_threads);
      
      % random Fourier coefficients
      fh = f_hat(rand((h.N+1)*(h.N+1),1));
      h.f_hat = fh;
      
      % Create plan of class NFSFT.
      plan = nfsft(N,M,bitor(NFSFT_NORMALIZED,NFSFT_USE_DPT));
      plan2 = nfsft(N,M,NFSFT_NORMALIZED);
      plan.x = X;
      plan2.x = X;

      plan.fhat = fh;
      plan2.fhat = fh;
          
      % NFSFT transform
      nfsft_trafo(plan);
      nfsft_trafo(plan2);
          
      % function values multiplied by quadrature weights
      h.f = plan.f .* W';
      h.f2 = plan2.f .* W';
    end
    
    function h = destroy(h)
      h.x = [];
      h.f_hat = [];
      h.f = [];
      h.f2 = [];
      h = destroy@nfsftTestcaseDelegate(h);
    end
  end
  
end

