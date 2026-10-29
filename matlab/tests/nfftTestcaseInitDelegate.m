classdef nfftTestcaseInitDelegate
  %NFFTTESTCASEDELEGATEINIT Summary of this class goes here
  %   Detailed explanation goes here

% Copyright (c) 2018, 2020 The NFFT Authors
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
    m = 0;
    flags = 0;
    fftw_flags = 0;
    K;
    sigma;
  end
  
  methods
    function h = nfftTestcaseInitDelegate(name, m, flags, fftw_flags)
      h.name = name;
      if (m==0)
        h.m = nfftmex('get_default_window_cut_off_m');
      else
        h.m = m;
      end
      h.flags = flags;
      h.fftw_flags = fftw_flags;
    end
    
    function [h, plan] = init(h, d, N, M)
      n = 2.^(ceil(log2(N))+1);
      h.sigma = n./N;
      switch h.name
        case 'init_guru'
          argsin = num2cell([d, N(:)', M, n(:)', h.m, h.flags, h.fftw_flags]);
          plan = nfft_init_guru(argsin{:});
        case 'init_class_no_flags'
          argsin = num2cell([n(:)', h.m]);
          plan = nfft(d,N,M,argsin{:});
        case 'init_class_flags'
          argsin = num2cell([n(:)', h.m, h.flags, h.fftw_flags]);
          plan = nfft(d,N,M,argsin{:});
        case 'init_1d'
          plan = nfft_init_1d(N, M);
        case 'init_2d'
          plan = nfft_init_2d(N(1), N(2), M);
        case 'init_3d'
          plan = nfft_init_3d(N(1), N(2), N(3), M);
        case 'init'
          plan = nfft_init(N, M);
        case 'init_class'
          plan = nfft(d, N, M);
        otherwise
          error('Unknown init name: %s', h.name);
      end
    end
  end
  
end

