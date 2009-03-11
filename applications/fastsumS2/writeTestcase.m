%
% Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
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

  % Write FPT threshold.
  fprintf(file,'threshold=%e\n',threshold);

end

% Write kernel type
fprintf(file,'kernel=%d\n',kernel);

% Write number of parameter sets.
fprintf(file,'parameter_sets=%d\n',size(parameters,1));

% Write number of parameters.
fprintf(file,'parameters=%d\n',size(parameters,2));

% Write parameters sets.
for j=1:size(parameters,1)
  for k=1:size(parameters,2)
    % Write parameter k of parameter sets j.
    fprintf(file,'%f\n',parameters(j,k));
  end
end

% Write number of bandwidths.
fprintf(file,'bandwidths=%d\n',length(bandwidths));

% Write bandwidths.
fprintf(file,'%d\n',bandwidths);

% Write number of node sets.
fprintf(file,'node_sets=%d\n',size(nodes,1));

% Write node sets.
for j = 1:size(nodes,1)
  fprintf(file,'L=%d\n',nodes(j,1));
  fprintf(file,'D=%d\n',nodes(j,2));
  fprintf(file,'compare=%d\n',nodes(j,3));
  if (nodes(j,3) == 1)
    fprintf(file,'precomputed=%d\n',nodes(j,4));
    fprintf(file,'repetitions=%d\n',nodes(j,5));
  end
end

% End of function
return;
