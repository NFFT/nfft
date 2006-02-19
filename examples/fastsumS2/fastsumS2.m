% fastsumS2 - Fast summation of radial functions on the sphere
infilename = 'data.in';
outfilename = 'testcase0.dat';
programname = 'fastsumS2';
selection = menu('fastsumS2 - Fast summation of radial functions on the sphere','5.1 (a)','5.1 (b)','5.1 (c)','5.1 (d)')
if (selection == 1)
  file = fopen(infilename,'w');
  % Number of testcases
  fprintf(file,'testcases=1\n');
  parameters = {};
  parameters{1,1,1} = 'h';
  parameters{1,1,2} = 0.8;
  parameters{2,1,1} = 'h';
  parameters{2,1,2} = 0.9;  nodes = {};
  nodes{1,1} = 10;
  nodes{1,2} = 10;
  nodes{1,3} = 1;
  nodes{1,4} = 0;
  nodes{1,5} = 1;
  nodes{2,1} = 100;
  nodes{2,2} = 100;
  nodes{2,3} = 1;
  nodes{2,4} = 0;
  nodes{2,5} = 1;
  parameters
  nodes
  writeTestcase(file,0,0,6,1000,0,parameters,4:4:64,nodes);
  fclose(file);
  system(sprintf('./%s < %s',programname,infilename));
  file = fopen(outfilename,'r');
  T = readTestcase(file)
%   kernel = fscanf(file,'kernel=%d\n')
%   nfsft = fscanf(file,'nfsft=%d\n')
%   if (nfsft == 1)
%     nfft = fscanf(file,'nfft=%d\n')
%     if (nfft == 1)
%       cutoff = fscanf(file,'cutoff=%d\n')
%     end
%     threshold = fscanf(file,'threshold=%d\n')
%   end
   fclose(file);
end
