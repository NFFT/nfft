% Copyright (c) 2018-2019 The NFFT Authors
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <https://www.gnu.org/licenses/>.

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
