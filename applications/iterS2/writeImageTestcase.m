function writeTestcase(file,usenfsft,usenfft,cutoff,usefpt,threshold,...
  bandwidth,im)

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

[m,n] = size(im);

r = length(find(im == 0))

theta = ((0:(m-1))+0.5)*pi/m;
phi = (0:(n-1))*2*pi/n;

% Write number of nodes.
fprintf(file,'nodes=%d\n',(m*n)-r);

% Write nodes and function values.
for j=1:length(theta)
  for k=1:length(phi)
    % Write node (j,k) and corresponding function value.
    if (abs(im(j,k)) ~= 0)
      fprintf(file,'%f %f %f %f\n',theta(j),phi(k),real(im(j,k)),imag(im(j,k)));
    end
  end
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
