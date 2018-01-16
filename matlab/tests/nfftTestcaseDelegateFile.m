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
      c = strfind(h.filename, '/');
      if ~isempty(c)
        c = c(end)+1;
      else
        c = 1;
      end
      fprintf('%-31s', h.filename(c:end));
      
      % Dimensions.
      h.d = fscanf(file, '%d', 1);
      % Bandwidths.
      h.N = zeros(h.d,1);
      for j=1:h.d
        h.N(j) = fscanf(file, '%d', 1);
      end
      % Number of nodes.
      h.M = fscanf(file, '%d', 1);
      
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
          x(j,t) = fscanf(file, '%g', 1);
        end
      end
      h.x = x;
      
      % Fourier coefficients.
      f_hat = zeros(h.NN,1);
      for j=1:h.NN
        val_real = fscanf(file, '%g', 1);
        val_imag = fscanf(file, '%g', 1);
        f_hat(j) = val_real + 1i*val_imag;
      end
      h.f_hat = f_hat;
      
      % Reference function values.
      f = zeros(h.M,1);
      for j=1:h.M
        val_real = fscanf(file, '%g', 1);
        val_imag = fscanf(file, '%g', 1);
        f(j) = val_real + 1i*val_imag;
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

