function writeTestcase(file,usenfsft,usenfft,cutoff,usefpt,threshold,...
  bandwidth,theta,phi,f)

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
