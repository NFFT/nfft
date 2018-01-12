%SUBSASGN Index assignment function for f_hat class
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
function p = subsasgn(p,s,v)
switch s.type
case '()'
  if (length(s.subs) > 2)
    error('Index must consist of two elements.');
  end

  ind1 = s.subs{1};
  ind2 = s.subs{2};

  if (any(round(ind1) ~= ind1) || min(ind1) < 0 || max(ind1) > p.N)
    error('Invalid index in first component');
  end

  if (any(round(ind2) ~= ind2) || min(ind2) < -p.N || max(ind2) > p.N)
    error('Invalid index in second component');
  end

  if (length(ind1) == 1 && min(ind2) >= -ind1 && max(ind2) <= ind1)
    p.f_hat(ind1^2+ind1+1+ind2) = v;
  elseif (length(ind2) == 1 && min(ind1) >= abs(ind2))
    ind = zeros(size(ind1));
    for k = 1:length(ind1)
      ind(k) = ind1(k)^2+ind1(k)+1+ind2;
    end
    p.f_hat(ind) = v;
  else
    error('Invalid subindex');
  end
otherwise
  error('Wrong subindex format');
end
