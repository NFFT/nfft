classdef nfsoftTestcaseTrafoDelegate
  %nfsoftTESTCASETRAFODELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    name = '';
  end
  
  methods
    function h = nfsoftTestcaseTrafoDelegate(name)
      h.name = name;
    end

    function plan = trafo(h, plan)
      switch h.name
        case 'trafo'
          if isnumeric(plan)
            nfsoft_trafo(plan);
          else
            plan.nfsoft_trafo;
          end
        case 'adjoint'
         if isnumeric(plan)
            nfsoft_adjoint(plan);
          else
            plan.nfsoft_adjoint;
          end
        otherwise
          error('Unknown trafo: %s', h.name);
      end
    end

    function result = isCost(h)
      switch h.name
%         case {'trafo_direct', 'adjoint_direct'}
%           result = 1;
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
    
    function val = acc(h, testcase)
      if testcase.N>20
        val = 2^30 * eps;
      else
        val = 2^22 * eps;
      end
    end
  end
  
end

