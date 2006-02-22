function C = readTestcase(file)
  tc_max = fscanf(file,'%d',1);
  C = cell(tc_max,1);
  for i=1:tc_max
    T = struct('usenfsft',[0],'usenfft',[0],'cutoff',[0],'usefpt',[0],...
      'threshold',[0],'kernel',[0],'parameters',{0},'bandwidths',[0],...
      'nodes',{0},'data',{0});
    % Use NFSFT
    v = fscanf(file,'%d',1);
    T.usenfsft = [v];
    if (T.usenfsft == 1)
      % Use NFFT
      v = fscanf(file,'%d',1);
      T.usenfft = [v];
      if (T.usenfft == 1)
        % NFFT cut-off parameter
        v = fscanf(file,'%d',1);
        T.cutoff = [v];
      end
      % Use FPT
      v = fscanf(file,'%d',1);
      T.usefpt = [v];
      % NFSFT threshold
      v = fscanf(file,'%lf',1);
      T.threshold = [v];
    end
    % Kernel type
    v = fscanf(file,'%d',1);
    T.kernel = [v];
    % Parameter sets
    ip_max = fscanf(file,'%d',1);
    ipp_max = fscanf(file,'%d',1);
    p = zeros(ipp_max,ip_max);
    p = fscanf(file,'%lf',[ipp_max,ip_max]);
    T.parameters = p';
    % Number of bandwidths
    bandwidths = fscanf(file,'%d',1);
    % Bandwidths
    m = fscanf(file,'%d',bandwidths);
    T.bandwidths = m;
    % Node sets
    ild_max = fscanf(file,'%d',1);
    nodes = zeros(ild_max,5);
    for j=1:ild_max
      % Node set
      nodes(j,1) = fscanf(file,'%d',1);
      nodes(j,2) = fscanf(file,'%d',1);
      nodes(j,3) = fscanf(file,'%d',1);
      if (nodes(j,3) == 1)
        nodes(j,4) = fscanf(file,'%d',1);
        nodes(j,5) = fscanf(file,'%d',1);
      end
    end
    datacell = {};
    for j=1:ip_max
      for k=1:ild_max
        data = fscanf(file,'%e',[6,bandwidths]);
        datacell{j,k} = data';
      end
    end
    T.data = datacell;
    T.nodes = nodes;
    C{i} = T;
  end

