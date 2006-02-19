% fastsumS2 - Fast summation of radial functions on the sphere
infilename = 'data.in';
outfilename = 'testcase0.dat';
programname = 'fastsumS2';
selection = menu('fastsumS2 - Fast summation of radial functions on the sphere','5.1 (a)','5.1 (b)','5.1 (c)','5.1 (d)')
if (selection == 1)
  file = fopen(infilename,'w');
  % Number of testcases
  fprintf(file,'testcases=1\n');
  % Use NFSFT
  fprintf(file,'nfsft=0\n');
  % Use NFFT
  %fprintf(file,'nfft=1\n');
  % NFFT cut-off parameter
  %fprintf(file,'cutoff=6\n');
  % NFSFT threshold
  %fprintf(file,'threshold=1e3\n');
  % Kernel type
  fprintf(file,'kernel=0\n');
  % Parameter sets
  fprintf(file,'parameter_sets=1\n');
  % Parameter h
  fprintf(file,'h=0.8\n');
  m = 4:4:64;
  % Number of bandwidths
  fprintf(file,'bandwidths=%d\n',length(m));
  % Bandwidths
  fprintf(file,'M=%d\n',m);
  % Node sets
  fprintf(file,'node_sets=1\n');
  fprintf(file,'L=1000\n');
  fprintf(file,'D=1000\n');
  fprintf(file,'compare=1\n');
  fprintf(file,'precomputed=0\n');
  fprintf(file,'repetitions=1\n');
  fclose(file);
  system(sprintf('./%s < %s',programname,infilename));
%   file = fopen(outfilename,'r');
%   kernel = fscanf(file,'kernel=%d\n')
%   nfsft = fscanf(file,'nfsft=%d\n')
%   if (nfsft == 1)
%     nfft = fscanf(file,'nfft=%d\n')
%     if (nfft == 1)
%       cutoff = fscanf(file,'cutoff=%d\n')
%     end
%     threshold = fscanf(file,'threshold=%d\n')
%   end
%   fclose(file);
end
  
  
  
  
  
  
