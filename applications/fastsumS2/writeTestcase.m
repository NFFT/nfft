function writeTestcase(file,usenfsft,usenfft,cutoff,usefpt,threshold,kernel,...
  parameters,bandwidths,nodes)

% WRITETESTCASE Write a fastsumS2 testcase definition to a file
%    WRITETESTCASE(FILE, USENFSFT, USENFFT, CUTOFF, USEFPT, THRESHOLD, KERNEL,
%    PARAMETERS, BANDWIDTHS, NODES) writes a fastsumS2 testcase specification to
%    the file associated with the file handle FILE.
%    The parameters are
%    - FILE The file handle associated with the file to write to,
%    - USENFSFT If true, the NFSFT algorithm is used, otherwise the direct NDSFT
%      algorithm,
%    - USENFFT If true, the NFFT algorithm is used, otherwise the direct NDSFT
%      algorithm (ignored if USENFSFT is false),
%    - CUTOFF The NFFT cut-off parameter (ignored if USENFFT is false)
%    - USEFPT If true, the fast polynomial transform is used, otherwise the
%      direct polynomial transform algorithm (ignored if USENFSFT is false),
%    - THRESHOLD The fast polynomial transform threshold parameter (ignored if
%      USEFPT is false),
%    - KERNEL The kernel function to be used (0 = Abel-Poisson kernel, 1 =
%      singularity kernel, 2 = locally supported kernel, 3 = spherical Gaussian
%      kernel),
%   - PARAMETERS A m x n matrix containing the kernel function parameters,
%     m is the number of parameter sets and n is the number of kernel function
%     parameters (1 for Abel-Poisson, singularity and spherical Gaussian kernel,
%     2 for the locally supported kernel)
%   - BANDWIDTH A vector containing the cut-off degrees for the approximation
%   - NODES A m x 5 matrix containing the node sets used, where m is the number
%     of different node sets. Each row contains a node set specification
%     containing
%     - in the first column the number of source nodes,
%     - in the second colum the number of target nodes,
%     - in the third column whether the direct sum evaluation has been performed
%       to compute the error E_infty,
%     - in the fourth column whether the precomputed direct sum evaluation has
%       been used (undefined if direct sum evaluation has not been used),
%     - in the fifth column the error E_infty (undefined if direct sum
%       evaluation has not been used).

% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%
% Write NFSFT usage flag.
fprintf(file,'nfsft=%d\n',usenfsft);

if (usenfsft == true)

  % Write NFFT usage flag.
  fprintf(file,'nfft=%d\n',usenfft);

  if (usenfft == true)

    % Write NFFT cut-off parameter.
    fprintf(file,'cutoff=%d\n',cutoff);

  end

  % Write FPT usage flag.
  fprintf(file,'fpt=%d\n',usefpt);

  % Write FPT threshold.
  fprintf(file,'threshold=%e\n',threshold);

end

% Write kernel type
fprintf(file,'kernel=%d\n',kernel);

% Write number of parameter sets.
fprintf(file,'parameter_sets=%d\n',size(parameters,1));

% Write number of parameters.
fprintf(file,'parameters=%d\n',size(parameters,2));

% Write parameters sets.
for j=1:size(parameters,1)
  for k=1:size(parameters,2)
    % Write parameter k of parameter sets j.
    fprintf(file,'%f\n',parameters(j,k));
  end
end

% Write number of bandwidths.
fprintf(file,'bandwidths=%d\n',length(bandwidths));

% Write bandwidths.
fprintf(file,'%d\n',bandwidths);

% Write number of node sets.
fprintf(file,'node_sets=%d\n',size(nodes,1));

% Write node sets.
for j = 1:size(nodes,1)
  fprintf(file,'L=%d\n',nodes(j,1));
  fprintf(file,'D=%d\n',nodes(j,2));
  fprintf(file,'compare=%d\n',nodes(j,3));
  if (nodes(j,3) == 1)
    fprintf(file,'precomputed=%d\n',nodes(j,4));
    fprintf(file,'repetitions=%d\n',nodes(j,5));
  end
end

% End of function
return;
