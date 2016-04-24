function writeTestcase(file,usenfsft,usenfft,cutoff,usefpt,threshold,...
  bandwidth,img,itheta,iphi)
%WRITEIMAGETESTCASE
%    WRITETESTCASE(FILE, USENFSFT, USENFFT, CUTOFF, USEFPT, THRESHOLD,
%    BANDWIDTH, THETA, PHI, F)
%
%   Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
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

% Generate Image
[m,n] = size(img)
dtheta = pi/m;
dphi = 2*pi/n;

if (exist('itheta','var') == false)
  ts = 1;
  te = m;
else
  ts = itheta(1)
  te = itheta(2)
end

if (exist('iphi','var') == false)
  ps = 1;
  pe = n;
else
  ps = iphi(1)
  pe = iphi(2)
end

im = img(ts:te,ps:pe);

r = length(find(im == 0))

theta = ((ts-0.5):(te-0.5))*dtheta;
phi = ((ps-1):(pe-1))*dphi;

% Write number of nodes.
fprintf(file,'nodes=%d\n',((te-ts+1)*(pe-ps+1))-r);

% Write nodes and function values.
for j=1:length(theta)
  for k=1:length(phi)
    % Write node (j,k) and corresponding function value.
    if (abs(im(j,k)) ~= 0)
      fprintf(file,'%f %f %f %f\n',theta(j),phi(k),real(im(j,k)),imag(im(j,k)));
    end
  end
end

theta = ((ts-0.5):0.25:(te-0.5))*dtheta;
phi = ((ps-1):0.25:(pe-1))*dphi;

% Write number of nodes.
fprintf(file,'nodes_eval=%d\n',length(theta)*length(phi));
% Write nodes and function values.
for j=1:length(theta)
  for k=1:length(phi)
    fprintf(file,'%f %f\n',theta(j),phi(k));
  end
end

% End of function
return;
