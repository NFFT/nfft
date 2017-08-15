function writeTestcase(file,usenfsft,usenfft,cutoff,usefpt,threshold,...
  testmode,gridtype,testfunction,repetitions,mode,bandwidths)
%WRITETESTCASE - Write qudratureS2 testcases
%
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
%
% Write usenfsft
fprintf(file,'nfsft=%d\n',usenfsft);
if (usenfsft == 1)
  % Write usenfft
  fprintf(file,'nfft=%d\n',usenfft);
  if (usenfft == 1)
    % Write NFFT cut-off parameter
    fprintf(file,'cutoff=%d\n',cutoff);
  end
  % Write use FPT
  fprintf(file,'fpt=%d\n',usefpt);
  if (usefpt == 1)
    % Write NFSFT threshold
    fprintf(file,'threshold=%e\n',threshold);
  end
end

fprintf(file,'testmode=%d\n',testmode);

if (testmode == 0)
  % Write grid type
  fprintf(file,'gridtype=%d\n',gridtype);
  % Write grid type
  fprintf(file,'testfunction=%d\n',testfunction(1));
  if (testfunction(1) == 0)
    fprintf(file,'bandlimit=%d\n',testfunction(2));
  end
  % Write number of repetitions
  fprintf(file,'repetitions=%d\n',repetitions);
  % Write mode
  fprintf(file,'mode=%d\n',mode(1));
  if (mode(1) == 1)
    % Write points
    fprintf(file,'points=%d\n',mode(2));
  end
end

% Write number of bandwidths
fprintf(file,'bandwidths=%d\n',size(bandwidths,2));

if (testmode == 0)
  % Write bandwidths
  fprintf(file,'%d %d\n',bandwidths);
else
  % Write bandwidths
  fprintf(file,'%d %d %d\n',bandwidths);
end

% Check if we need to provide also quadrature weights. This is the case if
% the Gauss-Legendre quadrature grid is used.
if (gridtype==0)
  writeWeights(file,bandwidths(2,:));
end
