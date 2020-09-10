clear all;
addpath ../nfft
ok = 1;

try
  tests = nfftUnitTests;
  
  result = nfft_check_1d_direct_file(tests); ok = min(ok, result);
  result = nfft_check_1d_fast_file(tests); ok = min(ok, result);
  result = nfft_check_adjoint_1d_direct_file(tests); ok = min(ok, result);
  result = nfft_check_adjoint_1d_fast_file(tests); ok = min(ok, result);
  result = nfft_check_1d_online(tests); ok = min(ok, result);
  result = nfft_check_adjoint_1d_online(tests); ok = min(ok, result);
  
  result = nfft_check_2d_direct_file(tests); ok = min(ok, result);
  result = nfft_check_2d_fast_file(tests); ok = min(ok, result);
  result = nfft_check_adjoint_2d_direct_file(tests); ok = min(ok, result);
  result = nfft_check_adjoint_2d_fast_file(tests); ok = min(ok, result);
  result = nfft_check_2d_online(tests); ok = min(ok, result);
  result = nfft_check_adjoint_2d_online(tests); ok = min(ok, result);
  
  result = nfft_check_3d_direct_file(tests); ok = min(ok, result);
  result = nfft_check_3d_fast_file(tests); ok = min(ok, result);
  result = nfft_check_adjoint_3d_direct_file(tests); ok = min(ok, result);
  result = nfft_check_adjoint_3d_fast_file(tests); ok = min(ok, result);
  result = nfft_check_3d_online(tests); ok = min(ok, result);
  result = nfft_check_adjoint_3d_online(tests); ok = min(ok, result);

  result = nfft_check_4d_online(tests); ok = min(ok, result);
  result = nfft_check_adjoint_4d_online(tests); ok = min(ok, result);
  
  % 5d tests only in double precision
  if (nfftmex('get_epsilon') < 1e-10)
    result = nfft_check_5d_online(tests); ok = min(ok, result);
    result = nfft_check_adjoint_5d_online(tests); ok = min(ok, result);
  end

  clear tests;
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
