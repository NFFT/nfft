classdef nfsftTestcaseInitDelegate
  %NFSFTTESTCASEDELEGATEINIT Summary of this class goes here
  %   Detailed explanation goes here

% Copyright (c) 2018 The NFFT Authors
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
  
  properties
    name = '';
    txt = '';
    flags = 0;
    kappa = 1000;
    fpt_flags = 0;
    nfft_flags = [];
    nfft_cutoff_m = [];
  end
  
  methods
    function h = nfsftTestcaseInitDelegate(name, flags, fpt_flags, nfft_flags, nfft_cutoff_m)
      h.name = name;
      h.txt = h.name;
      if exist('flags','var') && ~isempty(flags)
        h.flags = flags;
        if bitand(flags,NFSFT_USE_DPT) ~= 0
          h.txt = sprintf('%s DPT', h.name);
        end
      end
      if exist('fpt_flags','var') && ~isempty(fpt_flags)
        h.fpt_flags = fpt_flags;
      end
      if exist('nfft_flags','var') && ~isempty(nfft_flags)
        h.nfft_flags = nfft_flags;
      end
      if exist('nfft_cutoff_m','var') && ~isempty(nfft_cutoff_m)
        h.nfft_cutoff_m = nfft_cutoff_m;
      end
    end

    function nfsft_precompute(h, N)
      if ~isempty(h.fpt_flags)
        nfsft_precompute(N,h.kappa,h.flags,h.fpt_flags);
      elseif ~isempty(h.flags)
        nfsft_precompute(N,h.kappa,h.flags);
      elseif ~isempty(h.kappa)
        nfsft_precompute(N,h.kappa);
      else
        nfsft_precompute(N);
      end
    end
    
    function [h, plan] = init(h, N, M)
      switch h.name
        case 'init'
          h.nfsft_precompute(N);
          plan = nfsft_init(N,M);
        case 'init_advanced'
          h.nfsft_precompute(N);
          plan = nfsft_init_advanced(N,M,h.flags);
        case 'init_guru'
          h.nfsft_precompute(N);
          plan = nfsft_init_guru(N,M,h.flags,h.nfft_flags,h.nfft_cutoff_m);
        case 'init_class'
          if isempty(h.flags) && isempty(h.nfft_cutoff_m) && isempty(h.fpt_flags) && isempty(h.nfft_flags)
            plan = nfsft(N,M);
          elseif ~isempty(h.flags)
            plan = nfsft(N,M,h.flags);
          else
            plan = nfsft(N,M,h.flags,h.kappa,h.nfft_cutoff_m,h.fpt_flags,h.nfft_flags);
          end
        otherwise
          error('Unknown init name: %s', h.name);
      end
    end
  end
  
end

