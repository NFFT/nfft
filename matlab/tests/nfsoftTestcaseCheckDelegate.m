classdef nfsoftTestcaseCheckDelegate
  %UNTITLED3 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    name = '';
  end
  
  methods
    function h = nfsoftTestcaseCheckDelegate(name)
      h.name = name;
    end

    function plan = prepare(h, plan, f, f_hat)
      switch h.name
        case 'trafo'
          if isnumeric(plan)
            nfsoft_set_f_hat(plan, double(f_hat));
          else
            plan.fhat = f_hat;
          end
        case 'adjoint'
          if isnumeric(plan)
            nfsoft_set_f(plan, f);
          else
            plan.f = f;
          end
        otherwise
          error('Unknown check: %s', h.name);
      end
    end

    function result = compare(h, plan, f, f_hat)
      if isnumeric(plan)
        plan_f = nfsoft_get_f(plan);
        plan_f_hat = nfsoft_get_f_hat(plan);
      else
        plan_f = plan.f;
        plan_f_hat = plan.fhat;
      end
      switch h.name
        case 'trafo'
          numerator = max(abs(f - plan_f));
          denominator = sum(sum(abs(double(plan_f_hat))));
        case 'adjoint'
          numerator = max(max(abs(double(f_hat) - double(plan_f_hat))));
          denominator = sum(abs(plan_f));
        otherwise
          error('Unknown check: %s', h.name);
      end
      
      if numerator == 0
        result = 0;
      else
        result = numerator / denominator;
      end
    end
    
  end
  
end

