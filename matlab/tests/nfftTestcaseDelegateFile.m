classdef nfftTestcaseDelegateFile < nfftTestcaseDelegate
  %NFFTTESTCASEDELEGATE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Hidden=true,SetAccess='protected',GetAccess='public');
    filename = [];
  end
  
  methods
    function h = nfftTestcaseDelegateFile(filename)
%      h = h@nfftTestcaseDelegate(setup, destroy);
      h.filename = filename;
    end
    
    function h = setup(h)
      file = fopen(h.filename, 'r');
      if file == -1
        error('Cannot open file: %s', h.filename);
      end
      c = strfind(h.filename, '/');
      if ~isempty(c)
        c = c(end)+1;
      else
        c = 1;
      end
      fprintf('%-31s', h.filename(c:end));
      
      % Dimensions.
      [d,count] = fscanf(file, '%d', 1);
      if count ~= 1
        error('Unable to read d from file');
      end
      h.d = d;

      % Bandwidths.
      h.N = zeros(h.d,1);
      for j=1:h.d
        [N_j,count] = fscanf(file, '%d', 1);
        if count ~= 1
          error('Unable to read N from file');
        end
	h.N(j) = N_j;
      end

      % Number of nodes.
      [M,count] = fscanf(file, '%d', 1);
      if count ~= 1
        error('Unable to read M from file');
      end
      h.M = M;
      
      fprintf(' d = %-1d, N = [', h.d);
      if h.d > 3
        for j=1:h.d
          if j > 1
            fprintf(',');
          end
          fprintf('%-3d', h.N(j));
        end
        for j=1:4-h.d
          fprintf('%s%-3s', ' ', '');
        end
      else
        for j=1:h.d
          if j > 1
            fprintf(',');
          end
          fprintf('%-5d', h.N(j));
        end
        for j=1:3-h.d
          fprintf('%s%-5s', ' ', '');
        end
      end
      fprintf('],');
      fprintf(' M = %-5d', h.M);
      
      fprintf(' nthreads = %d', nfft_get_num_threads);

      h.NN = prod(h.N);
      
      % Nodes.
      x = zeros(h.M, h.d);
      for j=1:h.M
        for t=1:h.d
          [x_j_t,count] = fscanf(file, '%g', 1);
          if count ~= 1
            error('Unable to read x from file');
          end
          x(j,t) = x_j_t;
        end
      end
      h.x = x;
      
      % Fourier coefficients.
      f_hat = zeros(h.NN,1);
      for j=1:h.NN
        [val_real,count] = fscanf(file, '%g', 1);
        if count ~= 1
          error('Unable to read fhat from file');
        end
        [val_imag,count] = fscanf(file, '%g', 1);
        if count ~= 1
          error('Unable to read fhat from file');
        end
        f_hat(j) = val_real + 1i*val_imag;
      end
      h.f_hat = f_hat;
      
      % Reference function values.
      f = zeros(h.M,1);
      for j=1:h.M
        [val_real,count] = fscanf(file, '%g', 1);
        if count ~= 1
          error('Unable to read f from file');
        end
        [val_imag,count] = fscanf(file, '%g', 1);
        if count ~= 1
          error('Unable to read f from file');
        end
        f(j) = val_real + 1i*val_imag;
        if count ~= 1
          error('Unable to read fhat from file');
        end
      end
      h.f = f;
      
      fclose(file);
    end

    function h = destroy(h)
      h.d = [];
      h.N = [];
      h.NN = [];
      h.M = [];
      h.x = [];
      h.f_hat = [];
      h.f = [];
      h = destroy@nfftTestcaseDelegate(h);
    end
  end    
  
end

