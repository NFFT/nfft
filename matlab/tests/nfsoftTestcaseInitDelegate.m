classdef nfsoftTestcaseInitDelegate
  %nfsoftTESTCASEDELEGATEINIT Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    name = '';
    txt = '';
    flags = 0;
    kappa = 1000;
    nfft_flags = [];
    nfft_cutoff_m = [];
    fftw_size = [];
  end
  
  methods
    function h = nfsoftTestcaseInitDelegate(name, flags, nfft_flags, nfft_cutoff_m, kappa, fftw_size)
      h.name = name;
      h.txt = h.name;
      if exist('flags','var') && ~isempty(flags)
        h.flags = flags;
        if bitand(flags,NFSOFT_USE_DPT) ~= 0
          h.txt = sprintf('%s DPT', h.name);
        end
      end
      if exist('nfft_flags','var') && ~isempty(nfft_flags)
        h.nfft_flags = nfft_flags;
      end
      if exist('nfft_cutoff_m','var') && ~isempty(nfft_cutoff_m)
        h.nfft_cutoff_m = nfft_cutoff_m;
      end
      if exist('kappa','var') && ~isempty(kappa)
        h.kappa = kappa;
      end
      if exist('fftw_size','var') && ~isempty(fftw_size)
        h.fftw_size = fftw_size;
      end
    end
    
    function [h, plan] = init(h, N, M)
      if ~isempty(h.fftw_size)
        h.fftw_size = h.fftw_size(N);
      end
      switch h.name
        case 'init'
          if ~isempty(h.nfft_cutoff_m) || ~isempty(h.nfft_flags) || ~isempty(h.fftw_size)
            plan = nfsoft_init(N,M,h.flags,h.nfft_flags,h.nfft_cutoff_m,h.kappa,h.fftw_size);
          elseif ~isempty(h.flags)
            plan = nfsoft_init(N,M,h.flags);
          else
            plan = nfsoft_init(N,M);
          end
        case 'init_class'
          if ~isempty(h.nfft_cutoff_m) || ~isempty(h.nfft_flags) || ~isempty(h.fftw_size)
            plan = nfsoft(N,M,h.flags,h.nfft_flags,h.nfft_cutoff_m,h.kappa,h.fftw_size);
          elseif ~isempty(h.flags)
            plan = nfsoft(N,M,h.flags);
          else
            plan = nfsoft(N,M);
          end
        otherwise
          error('Unknown init name: %s', h.name);
      end
    end
  end
  
end

