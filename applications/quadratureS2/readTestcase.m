function T = readTestcase(file)
% READTESTCASE Read quadratureS2.c testcase results from file
%   READTESTCASE(FILE) reads fastsumS2.c testcase results from the file FILE.
%   The testcase results are returned as a cell vector containing the data for
%   each individual testcase as a structure with the following fields:

% $Id$

% Read the number of testcases.
tc_max = fscanf(file,'%d',1);

% Create empty cell array for the testcases.
T = cell(tc_max,1);

% Cycle through all testcases.
for i = 1:tc_max

  % Create structure.
  t = struct('USENFSFT',[0],'USENFFT',[0],'CUTOFF',[0],'USEFPT',[0],...
    'THRESHOLD',[0],'testmode',[0],'gridtype',{0},'testfunction',[0],...
    'bandlimit',[0],'repetitions',[0],'mode',{0},'points',[0],...
    'parameters',[0],'data',[0]);

  % Read NFSFT usage flag.
  v = fscanf(file,'%d',1);
  if (v >= 1)

    % Set NFSFT usage flag in the structure.
    t.USENFSFT = [true];

    % Read NFFT usage flag.
    v = fscanf(file,'%d',1);

    if (v >= 1)

      % Set NFSFT usage flag in the structure.
      t.USENFFFT = [true];

      % Read the NFFT cut-off parameter.
      v = fscanf(file,'%d',1);
      t.CUTOFF = [v];

    else

      % Set NFSFT usage flag in the structure.
      t.USENFFFT = [false];

    end

    % Read FPT usage flag.
    v = fscanf(file,'%d',1);

    if (v >= 1)

      % Set NFSFT usage flag in the structure.
      t.USEFPT = [true];

      % Read the FPT threshold.
      v = fscanf(file,'%lf',1);
      t.THRESHOLD = [v];
      
    else

      % Set NFSFT usage flag in the structure.
      t.USEFPT = [false];

    end

  else

    % Set NFSFT usage flag in the structure.
    t.USENFSFT = [false];

  end

  % Read the test mode.
  v = fscanf(file,'%d',1);
  t.testmode = [v];

  if (v == 0)

    % Read the grid type.
    v = fscanf(file,'%d',1);
    t.gridtype = [v];

    % Read the testfunction.
    v = fscanf(file,'%d',1);
    t.testfunction = [v];
    
    if (v == 0)

      % Read the badnlimit.
      v = fscanf(file,'%d',1);
      t.bandlimit = [v];

    end

    % Read the number of repetitions.
    v = fscanf(file,'%d',1);
    t.repetitions = [v];

    % Read the mode.
    v = fscanf(file,'%d',1);
    t.mode = [v];

    if (v == 1)

      % Read the points.
      v = fscanf(file,'%d',1);
      t.points = [v];

    end

  end
  
  % Read the number of bandwidths.
  v = fscanf(file,'%d',1);
  
  if (t.testmode == 0)
    % Create empty array for parameters.
    p = zeros(2,v);
    % Create empty array for parameters.
    data = zeros(3,v);

    % Read parameters.
    p = fscanf(file,'%d',[2,v]);

    % Transpose matrix to get dimensions right.
    t.parameters = p';

    % Read parameters.
    data = fscanf(file,'%lf',[3,v]);

    % Transpose matrix to get dimensions right.
    t.data = data';
  else
    % Create empty array for parameters.
    p = zeros(3,v);
    % Create empty array for parameters.
    data = zeros(1,v);

    % Read parameters.
    p = fscanf(file,'%d',[3,v]);

    % Transpose matrix to get dimensions right.
    t.parameters = p';

    % Read parameters.
    data = fscanf(file,'%lf',[1,v]);

    % Transpose matrix to get dimensions right.
    t.data = data';
  end

  % Assign testcase structure to field in cell array.
  T{i} = t;

end

% End of the function
return;
