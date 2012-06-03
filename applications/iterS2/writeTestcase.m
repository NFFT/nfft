% Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

% Write NFSFT usage flag.
fprintf(file,'nfsft=%d\n',usenfsft);

if (usenfsft == true)

  % Write NFFT usage flag.
  fprintf(file,'nfft=%d\n',usenfft);

  if (usenfft == true)

    % Write NFFT cut-off parameter.
    fprintf(file,'cutoff=%d\n',cutoff);

  end

  % Write FPT usage flag.
  fprintf(file,'fpt=%d\n',usefpt);

  if (usefpt == true)

    % Write FPT threshold.
    fprintf(file,'threshold=%e\n',threshold);

  end

end

% Write bandwidth
fprintf(file,'bandwidth=%d\n',bandwidth);

% Write number of nodes.
fprintf(file,'nodes=%d\n',length(theta));


% Write nodes and function values.
for j=1:length(theta)
  % Write node j and corresponding function value.
  fprintf(file,'%f %f %f %f\n',theta(j),phi(j),real(f(j)),imag(f(j)));
end

% Write number of nodes.
fprintf(file,'nodes_eval=%d\n',m*n);
% Write nodes and function values.
for j=1:length(theta)
  for k=1:length(phi)
    fprintf(file,'%f %f\n',theta(j),phi(k));
  end
end


% End of function
return;
