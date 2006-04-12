function writeTestcase(file,usenfsft,usenfft,cutoff,usefpt,threshold,kernel,parameters,bandwidths,nodes)
  % Use NFSFT
  fprintf(file,'nfsft=%d\n',usenfsft);
  if (usenfsft == 1)
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
  fprintf(file,'bandwidths=%d\n',length(bandwidths(:,1)));
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
