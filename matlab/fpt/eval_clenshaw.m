function [val] = eval_clenshaw(alpha, beta, gamma, a, nodes)
%EVAL_CLENSHAW Summary of this function goes here
%
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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
%

N = length(a)-1;

it2 = repmat(a(N+1),numel(nodes),1);
it1 = repmat(a(N),numel(nodes),1);
for k=N:-1:2
  temp = a(k-1) + it2 * gamma(k+1);
  it2 = it1 + it2 .* (alpha(k+1) * nodes + beta(k+1));
  it1 = temp;
end
it2 = it1 + it2 .* (alpha(2) * nodes + beta(2));
val = it2 * gamma(1);
end

