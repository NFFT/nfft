% fastsumS2 - Fast summation of radial functions on the sphere
infilename = 'data.in';
outfilename = 'testcase0.dat';
programname = 'fastsumS2';
selection = menu('fastsumS2 - Fast summation of radial functions on the sphere','Figure 5.1 (a)','Figure 5.1 (b)','Figure 5.1 (c)','Figure 5.1 (d)')
if (selection == 1)
  file = fopen(infilename,'w');
  % Number of testcases
  fprintf(file,'testcases=1\n');
  parameters = [0.8];
  nodes = [1000, 1000, 1, 0, 1];
  writeTestcase(file,1,1,6,1,1000,1,parameters,4:4:256,nodes);
  fclose(file);
  system(sprintf('./%s < %s',programname,infilename));
  file = fopen(outfilename,'r');
  T = readTestcase(file);
  createFigure(T.bandwidths,T.data{1,1}(:,6));
  fclose(file);
end
