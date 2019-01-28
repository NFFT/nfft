classdef nfsftTestcaseTrafoDelegate
  %NFSFTTESTCASETRAFODELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    name = '';
  end
  
  methods
    function h = nfsftTestcaseTrafoDelegate(name)
      h.name = name;
    end

    function plan = trafo(h, plan)
      switch h.name
        case 'trafo_direct'
          if isnumeric(plan)
            ndsft_trafo(plan);
          else
            plan.ndsft_trafo;
          end
        case 'trafo'
          if isnumeric(plan)
            nfsft_trafo(plan);
          else
            plan.nfsft_trafo;
          end
        case 'adjoint_direct'
          if isnumeric(plan)
            ndsft_adjoint(plan);
          else
            plan.ndsft_adjoint;
          end
        case 'adjoint'
         if isnumeric(plan)
            nfsft_adjoint(plan);
          else
            plan.nfsft_adjoint;
          end
        otherwise
          error('Unknown trafo: %s', h.name);
      end
    end

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
        otherwise
          error('not implemented');
      end
    end
    
    function val = acc(h, m)
      switch h.name
        case {'trafo_direct', 'adjoint_direct'}
          val = 2^13 * eps;
        otherwise
%           err = pi * (sqrt(m) + m) * sqrt(sqrt(1 - 1/2)) * exp(-2*pi * m * sqrt(1 - 1 / 2));
%           val = max(0.33 * err, 2500 * eps);
          val = 2^22 * eps;
      end
    end
  end
  
end

