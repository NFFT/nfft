function fastsumS2()
% $Id$

% The name of the input file
infilename = 'data.in';
% The name of the output file
outfilename = 'data.out';
% The name of the program
programname = 'quadratureS2';

% Display the menu.
selection = menu('quadratureS2 - Fast evaluation of quadrature formulae on the sphere',...
  'Generate Figure 6.1 (a)','Generate Figure 6.1 (b)','Generate Figure 6.2 (a)',...
  'Generate Figure 6.2 (b)')

% Open input data file.
file = fopen(infilename,'w');

if (selection == 1)
  % Set the grid type.
  % 0 = Gauss-Legendre
  gridtype = 0;
  % Set the number of repetitions.
  repetitions = 1;
  % Set the size parameters.
  S1 = 16:16:1024;
  %S1 = 16:16:128;
  S2 = [16,32,64,128,256,512,1024];
  %S2 = [16,32,64,128];
  % Set the bandwidhts.
  N1 = S1;
  N2 = S2;
  % Write the number of testcases.
  fprintf(file,'testcases=3\n');
  % Write the testcases.
  writeTestcase(file,1,1,3,1,1000,0,gridtype,[1],repetitions,[0],[N1;S1]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[1],repetitions,[0],[N1;S1]);
  writeTestcase(file,0,0,0,0,1000,0,gridtype,[1],repetitions,[0],[N2;S2]);
elseif (selection == 2)
  % Set the grid type.
  % 0 = Gauss-Legendre
  gridtype=0;
  % Set the number of repetitions.
  repetitions=1;
  % Set the bandwidhts.
  N=16:16:1024;
  % Set the size parameters.
  Q=1024*ones(size(N));
  % Write the number of testcases.
  fprintf(file,'testcases=4\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[3],repetitions,[0],[N;Q]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[4],repetitions,[0],[N;Q]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[5],repetitions,[0],[N;Q]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[6],repetitions,[0],[N;Q]);
elseif (selection == 3)
  % Set the grid type.
  % 0 = Gauss-Legendre
  gridtype=0;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  Q=16:16:1024;
  % Set the bandwidhts.
  N=128*ones(1,length(Q));
  % Write the number of testcases.
  fprintf(file,'testcases=1\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[0,128],repetitions,[0],[N;Q]);
elseif (selection == 4)
  % Set the grid type.
  % 1 = Clenshaw-Curtis
  gridtype=1;
  % Set the number of repetitions.
  repetitions=1;
  % Set the bandwidhts.
  N=16:16:1024;
  % Set the number of nodes.
  Q=N.*N;
  % Set the number of repetitions
  R=ceil(7500./(N.^(1.5))); %horzcat(100:-10:20,20*ones(1,length(N)-length(100:-10:20)));
  % Set the bandwidhts.
  N2=[1,2,4,8,16,32,64,128,256,512,1024];
  % Set the number of nodes.
  Q2=N2.*N2;
  % Set the number of repetitions
  R2=ceil(7500./(N2.^(2))); %horzcat(100:-10:20,20*ones(1,length(N)-length(100:-10:20)));
  % Write the number of testcases.
  fprintf(file,'testcases=2\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,1,gridtype,[0,128],repetitions,[0],[N;Q;R]);
  writeTestcase(file,0,1,6,1,1000,1,gridtype,[0,128],repetitions,[0],[N2;Q2;R2]);
elseif (selection == 5)
  % Set the grid type.
  % 0 = Gauss-Legendre
  gridtype=0;
  % Set the number of repetitions.
  repetitions=1;
  % Set the bandwidhts.
  N=720:16:1024;
  % Set the size parameters.
  Q=1024*ones(size(N));
  % Write the number of testcases.
  fprintf(file,'testcases=1\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[6],repetitions,[0],[N;Q]);
elseif (selection == 6)
  % Set the grid type.
  % 0 = Gauss-Legendre
  gridtype=0;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  %Q1=16:16:1024;
  %Q2=16:16:1024;
  Q3=[1024];
  % Set the bandwidhts.
  %N1=16:16:1024;
  %N2=16:16:1024;
  N3=[1024];
  % Write the number of testcases.
  fprintf(file,'testcases=1\n');
  % Write the testcases.
  %writeTestcase(file,1,1,3,1,1000,0,gridtype,[1],repetitions,[0],[N1;Q1]);
  %writeTestcase(file,1,1,6,1,1000,0,gridtype,[1],repetitions,[0],[N2;Q2]);
  writeTestcase(file,0,0,0,0,1000,0,gridtype,[1],repetitions,[0],[N3;Q3]);
elseif (selection == 7)
  % Set the grid type.
  % 3 = Equidiatribution
  gridtype=3;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  Q=16:16:2048;
  % Set the bandwidhts.
  N1=50*ones(1,length(Q));
  % Write the number of testcases.
  fprintf(file,'testcases=4\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,50],repetitions,[0],[N1;Q]);
elseif (selection == 8)
  % Set the grid type.
  % 3 = Equidiatribution
  gridtype=3;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  Q=16:16:2048;
  % Set the bandwidhts.
  N1=250*ones(1,length(Q));
  N2=500*ones(1,length(Q));
  N3=750*ones(1,length(Q));
  N4=1000*ones(1,length(Q));
  % Write the number of testcases.
  fprintf(file,'testcases=4\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,250],repetitions,[0],[N1;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,500],repetitions,[0],[N2;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,750],repetitions,[0],[N3;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,1024],repetitions,[0],[N4;Q]);
else
  error('Wrong selection!');
end

fclose(file);

if (selection == 1)
%if (false)
  system(sprintf('./%s < %s > %s',programname,infilename,outfilename));
  file = fopen(outfilename,'r');
  T = readTestcase(file);
  fclose(file);
  figure('Color',[1 1 1],'InvertHardcopy','off','PaperSize',[20.98 29.68]);
  axes('FontSize',16);
  x = T{1}.parameters(:,1);
  semilogy(x,T{1}.data(:,3),'-.','LineWidth',2,'Color',[0,0,0]);
  hold on
  semilogy(x,T{2}.data(:,3),'--','LineWidth',2,'Color',[0,0,0]);
  x2 = T{3}.parameters(:,1);
  semilogy(x2,T{3}.data(:,3),'-','LineWidth',2,'Color',[0,0,0]);
  semilogy(x2,T{3}.data(:,3),'.','MarkerSize',2,'Color',[0,0,0]);
  axis([x(1) x(end) 1e-16 1e-5])
  xlabel('S=M');
  ylabel('E_{\infty}','Rotation',0);
end
