classdef nfsoftUnitTests
  %UNITTEST Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    perform_exhaustive_tests_flag = 0;

    testcases_online_trafo = { ...
      nfsoftTestcaseDelegateOnline(4, 10, 'trafo') ...
      nfsoftTestcaseDelegateOnline(8, 10, 'trafo') ...
      nfsoftTestcaseDelegateOnline(16, 10, 'trafo') ...
                             };

    testcases_online_trafo_exhaustive = { ...
      nfsoftTestcaseDelegateOnline(16, 100, 'trafo') ...
      nfsoftTestcaseDelegateOnline(32, 10, 'trafo') ...
    };

    testcases_online_adjoint = { ...
      nfsoftTestcaseDelegateOnline(4, 10, 'adjoint') ...
      nfsoftTestcaseDelegateOnline(8, 10, 'adjoint') ...
      nfsoftTestcaseDelegateOnline(16, 10, 'adjoint') ...
                       };

    testcases_online_adjoint_exhaustive = { ...
      nfsoftTestcaseDelegateOnline(16, 100, 'adjoint') ...
      nfsoftTestcaseDelegateOnline(32, 10, 'adjoint') ...
    };

    initializers_std = { ...
%      nfsoftTestcaseInitDelegate('init') ...
%      nfsoftTestcaseInitDelegate('init',NFSOFT_USE_DPT,bitshift(1, 4),6) ...
      nfsoftTestcaseInitDelegate('init',NFSOFT_USE_DPT,[],[],[],@(N)3*N) ...
      nfsoftTestcaseInitDelegate('init',NFSOFT_USE_DPT,bitshift(1, 4),6,[],@(N)3*N) ...
      nfsoftTestcaseInitDelegate('init_class') ...
      nfsoftTestcaseInitDelegate('init_class', NFSOFT_USE_DPT) ...
    };

    testcases_gl_quadrature_online = { ...
      nfsoftTestcaseDelegateGLQuadratureOnline(4) ...
      nfsoftTestcaseDelegateGLQuadratureOnline(8) ...
      nfsoftTestcaseDelegateGLQuadratureOnline(16) ...
      nfsoftTestcaseDelegateGLQuadratureOnline(32) ...
    };

    testcases_gl_quadrature_online_exhaustive = { ...
      nfsoftTestcaseDelegateGLQuadratureOnline(64) ...
    };
  
    initializers_quadrature = { ...
      nfsoftTestcaseInitDelegate('init',NFSOFT_NORMALIZED,[],[],[],@(N)3*N) ...
      nfsoftTestcaseInitDelegate('init',bitor(NFSOFT_NORMALIZED,NFSOFT_USE_DPT),[],[],[],@(N)3*N) ...
      nfsoftTestcaseInitDelegate('init_class',NFSOFT_NORMALIZED) ...
      nfsoftTestcaseInitDelegate('init_class',bitor(NFSOFT_NORMALIZED,NFSOFT_USE_DPT)) ...
    };
  end
  
  methods
    function ok = check_single(~, testcase, init_delegate, check_delegate, trafo_delegate)
      ok = 0;
      
      testcase = testcase.setup;

      fprintf(', %-18s', init_delegate.txt);
      [init_delegate, plan] = init_delegate.init(testcase.N, testcase.M);

      if ~isempty(init_delegate.nfft_cutoff_m)
        fprintf(', m = %2d', init_delegate.nfft_cutoff_m);
      else
        fprintf(', m std.');
      end
      fprintf(', %-14s', trafo_delegate.name);
      if ~isempty(init_delegate.fftw_size)
        fprintf(', n = %3d', init_delegate.fftw_size);
      else
        fprintf(', n std. ');
      end

      if isnumeric(plan)
        nfsoft_set_x(plan, testcase.x);
      else
        plan.x = testcase.x;
      end

      run_test = 1;
      if (trafo_delegate.isCost)
        cost = trafo_delegate.cost(plan);
        if (cost > 0.1)
          fprintf(' -> %-4s (cost too high)\n','OK');
          ok = 1;
          run_test = 0;
        end
      end
      
      if run_test
        plan = check_delegate.prepare(plan, testcase.f, testcase.f_hat);
        
        plan = trafo_delegate.trafo(plan);
        
        err = check_delegate.compare(plan, testcase.f, testcase.f_hat);
        bound = trafo_delegate.acc(testcase);
        ok = err < bound;
        if ~ok; okstr = 'FAIL'; else, okstr = 'OK'; end
        fprintf(' -> %-4s %.3e (%.3e)\n', okstr, err, bound);
      end
      
      testcase.destroy;
      if isnumeric(plan)
        nfsoft_finalize(plan);
      else
        delete(plan);
      end
    end % function ok

    function result = check_many(h, testcases, initializers, check_delegate, trafos)
      result = 1;
      for k=1:length(trafos)
        for i=1:length(testcases)
          for j=1:length(initializers)
            r = h.check_single(testcases{i}, initializers{j}, check_delegate, trafos{k});
            result = min(result, r);
          end
        end
      end
    end

    function result = nfsoft_check_online(h)
      result = h.check_many(h.testcases_online_trafo, h.initializers_std, nfsoftTestcaseCheckDelegate('trafo'), {nfsoftTestcaseTrafoDelegate('trafo')});
      if h.perform_exhaustive_tests_flag
        result2 = h.check_many(h.testcases_online_trafo_exhaustive, h.initializers_std, nfsoftTestcaseCheckDelegate('trafo'), {nfsoftTestcaseTrafoDelegate('trafo')});
        result = min(result,result2);
      end
    end
    
    function result = nfsoft_check_adjoint_online(h)
      result = h.check_many(h.testcases_online_adjoint, h.initializers_std, nfsoftTestcaseCheckDelegate('adjoint'), {nfsoftTestcaseTrafoDelegate('adjoint')});
      if h.perform_exhaustive_tests_flag
        result2 = h.check_many(h.testcases_online_adjoint_exhaustive, h.initializers_std, nfsoftTestcaseCheckDelegate('adjoint'), {nfsoftTestcaseTrafoDelegate('adjoint')});
        result = min(result,result2);
      end
    end

    function ok = check_single_quadrature(~, testcase, init_delegate, trafo_delegate)
      ok = 0;
      
      init_txt = init_delegate.txt;
      fprintf('... %-18s', init_txt);
      [init_delegate, plan] = init_delegate.init(testcase.N, testcase.M);

      if ~isempty(init_delegate.nfft_cutoff_m)
        fprintf(', m = %2d', init_delegate.nfft_cutoff_m);
      else
        fprintf(', m std.');
      end
      if ~isempty(init_delegate.fftw_size)
        fprintf(', n = %3d', init_delegate.fftw_size);
      else
        fprintf(', n std. ');
      end
      fprintf(', %-14s', trafo_delegate.name);

      if isnumeric(plan)
        nfsoft_set_x(plan, testcase.x);
      else
        plan.x = testcase.x;
      end

      run_test = 1;
      if (trafo_delegate.isCost)
        cost = trafo_delegate.cost(plan);
        if (cost > 0.1)
          fprintf(' -> %-4s (cost too high)\n','OK');
          ok = 1;
          run_test = 0;
        end
      end
      
      if run_test
        fh = testcase.f_hat;

        f_mul_W = testcase.f;
        if isnumeric(plan)
          nfsoft_set_f(plan,f_mul_W);
        else
          plan.f = f_mul_W;
        end
          
        plan = trafo_delegate.trafo(plan);

        if isnumeric(plan)
          fh2 = nfsoft_get_f_hat(plan);
        else
          fh2 = plan.fhat;
        end

        err = norm(fh-fh2) / norm(fh);
        bound = 2^22* eps;
        ok = err < bound;
        if ~ok; okstr = 'FAIL'; else, okstr = 'OK'; end
        fprintf(' -> %-4s %.3e (%.3e)\n', okstr, err, bound);
      end
      
      if isnumeric(plan)
        nfsoft_finalize(plan);
      else
        delete(plan);
      end
    end
    
    function result = check_many_quadrature(h, testcases, initializers, trafos)
      result = 1;
      for i=1:length(testcases)
        testcase = testcases{i}.setup;
        bound = 2^22* eps;
        err = max(abs(testcase.f-testcase.f2)) / sum(sum(abs(double(testcase.f_hat))));
        ok = err < bound;
        if ~ok
          fprintf('Testcase setup FAILED\n');
          result = 0;
          continue;
        end

        fprintf('\n');

        for j=1:length(initializers)
          for k=1:length(trafos)
            r = h.check_single_quadrature(testcase, initializers{j}, trafos{k});
            result = min(result, r);
          end
        end
        testcase.destroy;
      end
    end

    function result = nfsoft_check_quadrature_online(h)
      result = h.check_many_quadrature(h.testcases_gl_quadrature_online, h.initializers_quadrature, {nfsoftTestcaseTrafoDelegate('adjoint')});
      if h.perform_exhaustive_tests_flag
        result2 = h.check_many_quadrature(h.testcases_gl_quadrature_online_exhaustive, h.initializers_quadrature, {nfsoftTestcaseTrafoDelegate('adjoint')});
        result = min(result,result2);
      end
    end

  end
end
