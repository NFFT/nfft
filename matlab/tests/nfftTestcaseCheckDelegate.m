classdef nfftTestcaseCheckDelegate
  %UNTITLED3 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    name = '';
  end
  
  methods
    function h = nfftTestcaseCheckDelegate(name)
      h.name = name;
    end

    function plan = prepare(h, plan, f, f_hat)
      switch h.name
        case 'trafo'
          if isnumeric(plan)
            nfft_set_f_hat(plan, f_hat);
          else
            if plan.d > 1
              N = plan.N(:)';
              f_hat = reshape(f_hat, N(end:-1:1));
              f_hat = permute(f_hat, plan.d:-1:1);
              f_hat = f_hat(:);
            end
            plan.fhat = f_hat;
          end
        case 'adjoint'
          if isnumeric(plan)
            nfft_set_f(plan, f);
          else
            plan.f = f;
          end
        otherwise
          error('Unknown check: %s', h.name);
      end
    end

    function result = compare(h, plan, f, f_hat)
      if isnumeric(plan)
        plan_f = nfft_get_f(plan);
        plan_f_hat = nfft_get_f_hat(plan);
      else
        plan_f = plan.f;
        plan_f_hat = plan.fhat;
       	if plan.d > 1
          N = plan.N(:)';
          plan_f_hat = reshape(plan_f_hat, N);
          plan_f_hat = permute(plan_f_hat, plan.d:-1:1);
          plan_f_hat = plan_f_hat(:);
	      end
      end
      switch h.name
        case 'trafo'
          numerator = max(abs(f - plan_f));
          denominator = sum(abs(plan_f_hat));
        case 'adjoint'
          numerator = max(abs(f_hat - plan_f_hat));
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

