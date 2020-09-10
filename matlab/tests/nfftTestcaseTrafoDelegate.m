classdef nfftTestcaseTrafoDelegate
  %UNTITLED3 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    name = '';
  end
  
  methods
    function h = nfftTestcaseTrafoDelegate(name)
      h.name = name;
    end

    function plan = trafo(h, plan)
      switch h.name
        case 'trafo_direct'
          if isnumeric(plan)
            ndft_trafo(plan);
          else
            plan.ndft_trafo;
          end
        case 'trafo'
          if isnumeric(plan)
            nfft_trafo(plan);
          else
            plan.nfft_trafo;
          end
        case 'adjoint_direct'
          if isnumeric(plan)
            ndft_adjoint(plan);
          else
            plan.ndft_adjoint;
          end
        case 'adjoint'
         if isnumeric(plan)
            nfft_adjoint(plan);
          else
            plan.nfft_adjoint;
          end
        otherwise
          error('Unknown trafo: %s', h.name);
      end
    end

%     function result = isCheck(h)
%       switch h.name
%         case {'trafo_direct', 'adjoint_direct'}
%           result = 0;
%         otherwise
%           result = 1;
%       end
%     end
%     
%     function ok = check(h, plan)
%       ok = 1;
% %       error('not implemented');
% %       switch h.name
% %         case {'trafo_direct', 'adjoint_direct'}
% %         otherwise
% %           %check
% %       end
%     end
    
    function result = isCost(h)
      switch h.name
        case {'trafo_direct', 'adjoint_direct'}
          result = 1;
        otherwise
          result = 0;
      end
    end
    
    function val = cost(h, plan)
      switch h.name
        case {'trafo_direct', 'adjoint_direct'}
          if isnumeric(plan)
            problem_size = length(nfft_get_f(plan)) * length(nfft_get_f_hat(plan));
          else
            problem_size = length(plan.f) * length(plan.fhat);
          end
          val = 1e-6 * problem_size;
        otherwise
          error('not implemented');
      end
    end
    
    function val = acc(h, m, sigma)
      epsilon = nfftmex('get_epsilon');
      switch h.name
        case {'trafo_direct', 'adjoint_direct'}
          val = 48 * epsilon;
        otherwise
          err = pi * (sqrt(m) + m) * sqrt(sqrt(1 - 1/2)) * exp(-2*pi * m * sqrt(1 - 1 / 2));
          val = max(0.33 * err, 2500 * epsilon);
      end
    end
  end
  
end

