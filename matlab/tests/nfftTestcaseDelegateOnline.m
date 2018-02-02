classdef nfftTestcaseDelegateOnline < nfftTestcaseDelegate
  %NFFTTESTCASEDELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Hidden=true,SetAccess='protected',GetAccess='public');
    trafo_type = '';
  end
  
  methods
    function h = nfftTestcaseDelegateOnline(d, N, M, trafo_type)
      h.d = d;
      h.N = N;
      h.M = M;
      h.NN = prod(N);
      h.trafo_type = trafo_type;
    end
    
    function h = setup(h)
      fprintf('%-31s', 'nfft_online');
      fprintf(' d = %-1d, N = [', h.d);
      if h.d > 3
        for j=1:h.d
          if j > 1
            fprintf(',');
          end
          fprintf('%-3d', h.N(j));
        end
        for j=1:4-h.d
          fprintf('%s%-3s', ' ', '');
        end
      else
        for j=1:h.d
          if j > 1
            fprintf(',');
          end
          fprintf('%-5d', h.N(j));
        end
        for j=1:3-h.d
          fprintf('%s%-5s', ' ', '');
        end
      end
      fprintf('],');
      fprintf(' M = %-5d', h.M);
      
      fprintf(' nthreads = %d', nfft_get_num_threads);
      
      % Nodes.
      h.x = rand(h.M, h.d) - 0.5;

      %
%       plan = nfft_init(h.N, h.M);
%       switch h.trafo_type
%         case 'trafo'
%           h.f_hat = rand(h.NN,1) - 0.5 + 1i * (rand(h.NN,1) - 0.5);
%           nfft_set_x(plan, h.x.');
%           nfft_set_f_hat(plan, h.f_hat);
%           ndft_trafo(plan);
%           h.f = nfft_get_f(plan);
%         case 'adjoint'
%           h.f = rand(h.M,1) - 0.5 + 1i * (rand(h.M,1) - 0.5);
%           nfft_set_x(plan, h.x.');
%           nfft_set_f(plan, h.f);
%           ndft_adjoint(plan);
%           h.f_hat = nfft_get_f_hat(plan);
%         otherwise
%           error('type %s not supported', h.trafo_type);
%       end
%       nfft_finalize(plan);
      if sum(mod(h.N,2)~=0) > 0
        error('only even N supported');
      end
      freq = (-h.N(h.d)/2:h.N(h.d)/2-1)';
      for t=h.d-1:-1:1
        freq = [reshape(repmat(-h.N(t)/2:h.N(t)/2-1, size(freq,1), 1), h.N(t)*size(freq,1), 1) repmat(freq, h.N(t), 1)];
      end
      switch h.trafo_type
        case 'trafo'
          h.f_hat = rand(h.NN,1) - 0.5 + 1i * (rand(h.NN,1) - 0.5);
          new_f = zeros(size(h.x,1),1);
          for j=1:size(h.x,1)
            new_f(j) = exp(-2*pi*1i*h.x(j,:)*freq') * h.f_hat;
          end
          h.f = new_f;
        case 'adjoint'
          h.f = rand(h.M,1) - 0.5 + 1i * (rand(h.M,1) - 0.5);
          new_f_hat = zeros(h.NN,1);
          for j=1:size(h.x,1)
            new_f_hat = new_f_hat + exp(2*pi*1i*freq*h.x(j,:)') * h.f(j);
          end
          h.f_hat = new_f_hat;
        otherwise
          error('type %s not supported', h.trafo_type);
      end

    end
    
    function h = destroy(h)
      h.x = [];
      h.f_hat = [];
      h.f = [];
      h = destroy@nfftTestcaseDelegate(h);
    end
  end
  
end

