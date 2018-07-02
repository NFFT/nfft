addpath ../nfsft
ok = 1;
bound = 2^20* eps;
bound_direct = 2^12* eps;
bound_direct_adjoint = 2^13* eps;
bound_reconstruction = 2^22* eps;

fprintf('Number of threads: %d\n', nfsft_get_num_threads());

try
  for N = 2.^(2:10)
    % precomputation
    nfsft_precompute(N,1000);

    % Comparison with accurate (slow) algorithm from Matlab
    % random Fourier coefficients
    fh = f_hat(rand((N+1)*(N+1),1)./(1:(N+1)*(N+1))');

    % number of nodes
    M = 80;

    % random nodes
    ph=rand(1,M)*2*pi;
    th=rand(1,M)*pi;
    X=[ph;th];

    % Create plan.
    plan = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT));

    % Set nodes.
    nfsft_set_x(plan,X);

    % set Fourier coefficients
    nfsft_set_f_hat(plan,double(fh));

    % direct (slow) transform
    ndsft_trafo(plan);
    f_direct = nfsft_get_f(plan);

    % fast NFSFT transform
    nfsft_trafo(plan);
    f_fast = nfsft_get_f(plan);


    % random function values for adjoint
    f = rand(M,1);

    % adjoint (slow) transform
    nfsft_set_f(plan,f)
    ndsft_adjoint(plan);
    fh2_direct = f_hat(nfsft_get_f_hat(plan));

    % Fast adjoint transform
    nfsft_set_f(plan,f)
    nfsft_adjoint(plan);
    fh2_fast = f_hat(nfsft_get_f_hat(plan));

    % finalize plan
    nfsft_finalize(plan);

    % the same with DPT
    plan_dpt = nfsft_init_advanced(N,M,bitor(NFSFT_NORMALIZED,NFSFT_PRESERVE_F_HAT)+NFSFT_USE_DPT);
    nfsft_set_x(plan_dpt,X);
    nfsft_set_f_hat(plan_dpt,double(fh));
    nfsft_trafo(plan_dpt);
    f_dpt = nfsft_get_f(plan_dpt);

    nfsft_set_f(plan_dpt,f)
    nfsft_adjoint(plan_dpt);
    fh2_dpt = f_hat(nfsft_get_f_hat(plan_dpt));
    nfsft_finalize(plan_dpt);


    % Algorithm in Matlab
    % this is very slow in Octave
    if N<200 || ~exist('octave_config_info')
      f_mat = legendre(0,cos(X(2,:)),'norm').' / sqrt(2*pi) * fh(0,0);
      fh2_mat = zeros(size(fh.f_hat));
      fh2_mat(1)= legendre(0,cos(X(2,:)),'norm') / sqrt(2*pi) * f;
      for n=1:N
        Lp = legendre(n,cos(X(2,:)),'norm') / sqrt(2*pi);
        f_mat=f_mat + (flipud(Lp).* exp((-n:0)'.*1i.*X(1,:))).'* fh(n,-n:0)...
         + (Lp(2:end,:).* exp((1:n)'.*1i.*X(1,:))).'* fh(n,1:n);

        fh2_mat(n^2+1:n^2+n+1)= (flipud(Lp).* exp(-(-n:0)'.*1i.*X(1,:)))* f;
        fh2_mat(n^2+n+2:n^2+2*n+1)= (Lp(2:end,:).* exp(-(1:n)'.*1i.*X(1,:)))* f;
      end
      fh2_mat = f_hat(fh2_mat);

      % Compute error
      error_trafo = [norm(f_fast-f_mat), norm(f_dpt-f_mat), norm(f_direct-f_mat)]/norm(f_mat);
      if (max(error_trafo) > bound) || (error_trafo(3) > bound_direct)
        ok = 0;
        result = 'FAIL';
      else
        result = ' OK ';
      end
      fprintf('NFSFT_trafo    (N =%4d) ---> %s  fpt: %.3e, dpt: %.3e (%.3e), direct: %.3e (%.3e)\n', ...
        N, result, error_trafo(1), error_trafo(2), bound, error_trafo(3), bound_direct);

      error_adjoint = [norm(fh2_fast-fh2_mat), norm(fh2_dpt-fh2_mat), norm(fh2_direct-fh2_mat)]/norm(fh2_mat);
      if (max(error_adjoint) > bound) || (error_adjoint(3) > bound_direct_adjoint)
        ok = 0;
        result = 'FAIL';
      else
        result = ' OK ';
      end
      fprintf('NFSFT_adjoint  (N =%4d) ---> %s  fpt: %.3e, dpt: %.3e (%.3e), direct: %.3e (%.3e)\n', ...
        N, result, error_adjoint(1), error_adjoint(2), bound, error_adjoint(3), bound_direct_adjoint);
    
    else %Octave for high N
      error_trafo = [norm(f_fast-f_direct), norm(f_dpt-f_direct)]/norm(f_direct);
      if (max(error_trafo) > bound)
        ok = 0;
        result = 'FAIL';
      else
        result = ' OK ';
      end
      fprintf('NFSFT_trafo    (N =%4d) ---> %s  fpt: %.3e, dpt: %.3e (%.3e)\n', ...
        N, result, error_trafo(1), error_trafo(2), bound);

      error_adjoint = [norm(fh2_fast-fh2_direct), norm(fh2_dpt-fh2_direct)]/norm(fh2_direct);
      if (max(error_adjoint) > bound) 
        ok = 0;
        result = 'FAIL';
      else
        result = ' OK ';
      end
      fprintf('NFSFT_adjoint  (N =%4d) ---> %s  fpt: %.3e, dpt: %.3e (%.3e)\n', ...
        N, result, error_adjoint(1), error_adjoint(2), bound);
    end

    % Run test with many nodes and Gauss-Legendre quadrature
    [X,W] = gl(N);
    M = size(X,2);

    for use_dpt = 1:2
      % Create plan of class NFSFT.
      plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED + (use_dpt-1)*NFSFT_USE_DPT);
      nfsft_set_x(plan,X);

      % random Fourier coefficients
      fh = f_hat(rand((N+1)*(N+1),1));
      nfsft_set_f_hat(plan,double(fh));

      % NFSFT transform
      nfsft_trafo(plan);

      % function values
      f = nfsft_get_f(plan);

      % adjoint transform, using quadrature weights to recover Fourier coefficients
      % the adjoint is only the inverse with the right quadrature weights W
      nfsft_set_f(plan,f.*W');
      nfsft_adjoint(plan);

      fh2 = nfsft_get_f_hat(plan);
      nfsft_finalize(plan)
      error_quad(use_dpt) = norm(fh2-fh)/norm(fh);
    end

    if max(error_quad) > bound_reconstruction
      ok=0;
      result='FAIL';
    else
      result=' OK ';
    end
    fprintf('Reconstruction (N =%4d) ---> %s  fpt: %.3e  dpt: %.3e (%.3e)\n',...
      N,result,error_quad(1),error_quad(2),bound_reconstruction);
  end
  
  nfsft_forget();

catch err
  try
    fprintf('Exception %s %s\n', err.identifier, err.message);
    err
  catch
  end
  ok = 0;
end

clear result;

if ok ~= 1
  fprintf('nfftUnitTest: at least one test failed\n');
  exit(1);
  return;
end
fprintf('nfftUnitTest: all tests succeeded\n');
exit(0);
