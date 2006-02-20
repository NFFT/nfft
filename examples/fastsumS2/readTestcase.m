function T = readTestcase(file)
  T = struct('usenfsft',[0],'usenfft',[0],'cutoff',[0],'threshold',[0],...
    'kernel',[0],'parameters',{0},'bandwidths',[0],'nodes',{0});
  v = 0;
  % Use NFSFT
  fscanf(file,'nfsft=%d\n',v);
  T.usenfsft = [v]; 
  if (T.usenfsft == 1)
    % Use NFFT
    fscanf(file,'nfft=%d\n',v);
    T.usenfft = [v]; 
    if (T.usenfft == 1)
      % NFFT cut-off parameter
      fscanf(file,'cutoff=%d\n',v);
      T.cutoff = [v];
    end  
    % NFSFT threshold
    fscanf(file,'threshold=%e\n',v);
    T.threshold = [v];
  end  
  % Kernel type
  fscanf(file,'kernel=%d\n',v);
  T.kernel = [v];
%  % Parameter sets
%  N = 0;
%  fscanf(file,'parameter_sets=%d\n',N);
%   for j=1:N
%     for k=1:size(parameters,2)
%       % Parameter h
%       fprintf(file,'%s=%f\n',parameters{j,k,1},parameters{j,k,2});
%     end
%   end  
%   % Number of bandwidths
%   fprintf(file,'bandwidths=%d\n',length(bandwidths));
%   % Bandwidths
%   fprintf(file,'M=%d\n',bandwidths);
%   % Node sets
%   fprintf(file,'node_sets=%d\n',size(nodes,1));
%   for j=1:size(nodes,1)
%     % Node set
%     fprintf(file,'L=%d\n',nodes{j,1});
%     fprintf(file,'D=%d\n',nodes{j,2});
%     fprintf(file,'compare=%d\n',nodes{j,3});
%     if (nodes{j,3} == 1)
%       fprintf(file,'precomputed=%d\n',nodes{j,4});
%       fprintf(file,'repetitions=%d\n',nodes{j,5});
%     end
%   end  
