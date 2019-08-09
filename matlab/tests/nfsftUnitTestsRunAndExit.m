clear all;
addpath ../nfsft
ok = 1;

fprintf('Number of threads: %d\n', nfsft_get_num_threads());

try
  tests = nfsftUnitTests;

  result = tests.nfsft_check_online; ok = min(ok, result);
  result = tests.nfsft_check_adjoint_online; ok = min(ok, result);

  result = tests.nfsft_check_trafo_equispaced_online; ok = min(ok, result);
  result = tests.nfsft_check_adjoint_equispaced_online; ok = min(ok, result);
  result = tests.nfsft_check_quadrature_online; ok = min(ok, result);
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
  fprintf('nfsftUnitTest: at least one test failed\n');
  exit(1);
  return;
end
fprintf('nfsftUnitTest: all tests succeeded\n');
exit(0);
