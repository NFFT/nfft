classdef nfsoftTestcaseDelegateOnline < nfsoftTestcaseDelegate
  %nfsoftTESTCASEDELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Hidden=true,SetAccess='protected',GetAccess='public')
    trafo_type = '';
  end
  
  methods
    function h = nfsoftTestcaseDelegateOnline(N, M, trafo_type)
      h.N = N;
      h.M = M;
      h.trafo_type = trafo_type;
    end
    
    function h = setup(h)
      fprintf('%-31s', 'nfsoft_online');
      fprintf(' N = %-5d,', h.N);
      fprintf(' M = %-5d,', h.M);
      
      fprintf(' nthreads = %d', nfsoft_get_num_threads);
      
      % random nodes
      alpha = rand(h.M,1)*2*pi;
      beta = acos(rand(h.M,1)*2-1);
      gamma = rand(h.M,1)*2*pi;
      X = [alpha';beta';gamma'];
      h.x = X;

     switch h.trafo_type
      case 'trafo'
        fh = rand(nfsoft_f_hat_size(h.N),1);
        h.f_hat = fh;
        f_mat = zeros(h.M,1);
        for m=1:h.M
            for n=0:h.N
                D = wignerD(n,alpha(m),beta(m),gamma(m));
                D = D(:);
                f_mat(m) = f_mat(m) + D.' * fh(nfsoft_f_hat_size(n-1)+1:nfsoft_f_hat_size(n));
            end
        end
        h.f = f_mat;
      case 'adjoint'
        f = rand(h.M,1) - 0.5 + 1i * (rand(h.M,1) - 0.5);
        h.f = f;
        fh2_mat = zeros(nfsoft_f_hat_size(h.N),1);
        for m=1:h.M
            for n=0:h.N
                D = wignerD(n,alpha(m),beta(m),gamma(m));
                D = D(:);
                fh2_mat(nfsoft_f_hat_size(n-1)+1:nfsoft_f_hat_size(n)) = ...
                  fh2_mat(nfsoft_f_hat_size(n-1)+1:nfsoft_f_hat_size(n)) + conj(D) * f(m);
            end
        end
        h.f_hat = fh2_mat;
        otherwise
          error('type %s not supported', h.trafo_type);
      end
    end
    
    function h = destroy(h)
      h.x = [];
      h.f_hat = [];
      h.f = [];
      h = destroy@nfsoftTestcaseDelegate(h);
    end
  end
  
end

