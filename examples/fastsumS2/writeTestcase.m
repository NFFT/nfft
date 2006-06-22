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
%    - PARAMETERS A vector containing the kernel parameters (1 for Abel-Poisson,
%      singularity and spherical Gaussian kernel, 2 for the locally supported
%      kernel),
%    - BANDWIDTHS A vector of cut-off degrees for the approximation
%    - NODES The nodes.

% Use NFSFT
fprintf(file,'nfsft=%d\n',usenfsft);
if (usenfsft >= 1)
  % Use NFFT
  fprintf(file,'nfft=%d\n',usenfft);
  if (usenfft == 1)
    % NFFT cut-off parameter
    fprintf(file,'cutoff=%d\n',cutoff);
  end
  % Use FPT
  fprintf(file,'fpt=%d\n',usefpt);
  % NFSFT threshold
  fprintf(file,'threshold=%e\n',threshold);
end
% Kernel type
fprintf(file,'kernel=%d\n',kernel);
% Parameter sets
fprintf(file,'parameter_sets=%d\n',size(parameters,1));
% Parameter sets
fprintf(file,'parameters=%d\n',size(parameters,2));
for j=1:size(parameters,1)
  for k=1:size(parameters,2)
    % Parameter h
    fprintf(file,'%f\n',parameters(j,k));
  end
end
% Number of bandwidths
bandwidths
fprintf(file,'bandwidths=%d\n',length(bandwidths));
% Bandwidths
fprintf(file,'%d\n',bandwidths);
% Node sets
fprintf(file,'node_sets=%d\n',size(nodes,1));
for j=1:size(nodes,1)
  % Node set
  fprintf(file,'L=%d\n',nodes(j,1));
  fprintf(file,'D=%d\n',nodes(j,2));
  fprintf(file,'compare=%d\n',nodes(j,3));
  if (nodes(j,3) == 1)
    fprintf(file,'precomputed=%d\n',nodes(j,4));
    fprintf(file,'repetitions=%d\n',nodes(j,5));
  end
end
