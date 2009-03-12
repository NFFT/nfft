%
% Get f_hat properties from the specified object
%
% Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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
% $Id$
function val = get(p, propName)
% GET Get f_hat properties from the specified object
% and return the value
%
%   Copyright (c) 2006, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
switch propName
case 'Data'
    val = p.f_hat;
otherwise
    error([propName,' Is not a valid f_hat property']);
end
