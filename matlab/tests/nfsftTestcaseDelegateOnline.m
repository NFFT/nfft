classdef nfsftTestcaseDelegateOnline < nfsftTestcaseDelegate
  %NFSFTTESTCASEDELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Hidden=true,SetAccess='protected',GetAccess='public')
    trafo_type = '';
  end
  
  methods
    function h = nfsftTestcaseDelegateOnline(N, M, trafo_type)
      h.N = N;
      h.M = M;
      h.trafo_type = trafo_type;
    end
    
    function h = setup(h)
      
      if (h.M==0)  % equispaced nodes
        fprintf('%-31s', 'nfsft_online_equispaced');
        ph=(-h.N-1:h.N)/(2*h.N+2)*2*pi;
        th=(0:h.N+1)/(2*h.N+2)*2*pi;
        [ph,th]=meshgrid(ph,th);
        h.x=[ph(:)';th(:)'];
        h.M=size(h.x,2);
      else % random nodes
        fprintf('%-31s', 'nfsft_online');
        ph=rand(1,h.M)*2*pi;
        th=rand(1,h.M)*pi;
        h.x=[ph;th];
      end
      
      fprintf(' N = %-5d,', h.N);
      fprintf(' M = %-5d,', h.M);
      
      fprintf(' nthreads = %d', nfsft_get_num_threads);

      switch h.trafo_type
        case 'trafo'
          fh = f_hat(rand((h.N+1)*(h.N+1),1));
          h.f_hat = fh;
          
          p = nfsft(h.N,h.M,NFSFT_USE_DPT,1000.0,6,FPT_NO_FAST_ALGORITHM);
          p.fhat = fh;
          p.x = h.x;
          nfsft_trafo_direct(p);
          h.f = p.f;
        case 'adjoint'
          f = rand(h.M,1) - 0.5 + 1i * (rand(h.M,1) - 0.5);
          h.f = f;

          p = nfsft(h.N,h.M,NFSFT_USE_DPT,1000.0,6,FPT_NO_FAST_ALGORITHM);
          p.f = f;
          p.x = h.x;
          nfsft_adjoint_direct(p);
          h.f_hat = p.fhat;
        otherwise
          error('type %s not supported', h.trafo_type);
      end

    end
    
    function h = destroy(h)
      h.x = [];
      h.f_hat = [];
      h.f = [];
      h = destroy@nfsftTestcaseDelegate(h);
    end
  end
  
end

