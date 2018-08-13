addpath ../nfsoft ../nfsft
ok = 1;

fprintf('Number of threads: %d\n', nfsoft_get_num_threads());

try
  tests = nfsoftUnitTests;
  if exist('perform_exhaustive_tests_flag','var')
    tests.perform_exhaustive_tests_flag = perform_exhaustive_tests_flag;
  end

  result = tests.nfsoft_check_online; ok = min(ok, result);
  result = tests.nfsoft_check_adjoint_online; ok = min(ok, result);

  result = tests.nfsoft_check_quadrature_online; ok = min(ok, result);
catch err
  try
    fprintf('Exception %s %s\n', err.identifier, err.message);
    disp(err)
  catch
  end
  ok = 0;
end

clear result;

if ok ~= 1
  fprintf('nfsoftUnitTest: at least one test failed\n');
  exit(1);
  return;
end
fprintf('nfsoftUnitTest: all tests succeeded\n');
exit(0);
