function writeTestcase(file,usenfsft,usenfft,cutoff,usefpt,threshold,...
  bandwidth,img,itheta,iphi)

% WRITETESTCASE Write a iterS2 testcase definition to a file
%    WRITETESTCASE(FILE, USENFSFT, USENFFT, CUTOFF, USEFPT, THRESHOLD,
%    BANDWIDTH, THETA, PHI, F)
%
% $Id$
%
% Copyright (C) 2005 Jens Keiner

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
