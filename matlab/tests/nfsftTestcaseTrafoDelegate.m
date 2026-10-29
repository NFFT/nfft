classdef nfsftTestcaseTrafoDelegate
  %NFSFTTESTCASETRAFODELEGATE Summary of this class goes here
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
  end
  
  methods
    function h = nfsftTestcaseTrafoDelegate(name)
      h.name = name;
    end

    function plan = trafo(h, plan)
      switch h.name
        case 'trafo_direct'
          if isnumeric(plan)
            ndsft_trafo(plan);
          else
            plan.ndsft_trafo;
          end
        case 'trafo'
          if isnumeric(plan)
            nfsft_trafo(plan);
          else
            plan.nfsft_trafo;
          end
        case 'adjoint_direct'
          if isnumeric(plan)
            ndsft_adjoint(plan);
          else
            plan.ndsft_adjoint;
          end
        case 'adjoint'
         if isnumeric(plan)
            nfsft_adjoint(plan);
          else
            plan.nfsft_adjoint;
          end
        otherwise
          error('Unknown trafo: %s', h.name);
      end
    end

    function result = isCost(h)
      switch h.name
        case {'trafo_direct', 'adjoint_direct'}
          result = 1;
        otherwise
          result = 0;
      end
    end
    
    function val = cost(h, plan)
      switch h.name
        otherwise
          error('not implemented');
      end
    end
    
    function val = acc(h, m)
      switch h.name
        case {'trafo_direct', 'adjoint_direct'}
          val = 2^13 * eps;
        otherwise
%           err = pi * (sqrt(m) + m) * sqrt(sqrt(1 - 1/2)) * exp(-2*pi * m * sqrt(1 - 1 / 2));
%           val = max(0.33 * err, 2500 * eps);
          val = 2^22 * eps;
      end
    end
  end
  
end

