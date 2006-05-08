function fastsumS2()
% $Id$
%
% quadratureS2 - Fast evaluation of quadrature formulae on the sphere
%
% Copyright (C) 2005 Jens Keiner
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

% The name of the input file
infilename = 'data.in';
% The name of the output file
outfilename = 'data.out';
% The name of the program
programname = 'quadratureS2';

% Display the menu.
selection = menu('quadratureS2 - Fast evaluation of quadrature formulae on the sphere',...
  'Accuracy','Truncation f3-f6','Approximation 128 GL',...
  'Approximation 128 CC','Approximation 128 HP','Approximation 128 ED','Equidistribution small',...
  'Equidistribution')

% Open input data file.
file = fopen(infilename,'w');

if (selection == 1)
  % Set the grid type.
  % 0 = Gauss-Legendre
  gridtype=0;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  %Q1=16:16:1024;
  %Q2=16:16:1024;
  Q3=[512];
  % Set the bandwidhts.
  %N1=16:16:1024;
  %N2=16:16:1024;
  N3=[512];
  % Write the number of testcases.
  fprintf(file,'testcases=1\n');
  % Write the testcases.
  %writeTestcase(file,1,1,3,1,1000,0,gridtype,[1],repetitions,[N1;Q1]);
  %writeTestcase(file,1,1,6,1,1000,0,gridtype,[1],repetitions,[N2;Q2]);
  writeTestcase(file,0,0,0,0,1000,0,gridtype,[1],repetitions,[N3;Q3]);
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
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[3],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[4],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[5],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[6],repetitions,[N;Q]);
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
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[0,128],repetitions,[N;Q]);
elseif (selection == 4)
  % Set the grid type.
  % 1 = Clenshaw-Curtis
  gridtype=1;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  Q=16:16:1024;
  % Set the bandwidhts.
  N=128*ones(1,length(Q));
  % Write the number of testcases.
  fprintf(file,'testcases=1\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[0,128],repetitions,[N;Q]);
elseif (selection == 5)
  % Set the grid type.
  % 2 = HEALPix
  gridtype=2;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  Q=[1,2,4,8,16,32,64,128,256,512];
  % Set the bandwidhts.
  N=128*ones(1,length(Q));
  % Write the number of testcases.
  fprintf(file,'testcases=1\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[0,128],repetitions,[N;Q]);
elseif (selection == 6)
  % Set the grid type.
  % 4 = Equidistribution
  gridtype=4;
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  Q=16:16:1024;
  % Set the bandwidhts.
  N=128*ones(1,length(Q));
  % Write the number of testcases.
  fprintf(file,'testcases=1\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[0,128],repetitions,[N;Q]);
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
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,50],repetitions,[N1;Q]);
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
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,250],repetitions,[N1;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,500],repetitions,[N2;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,750],repetitions,[N3;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,1024],repetitions,[N4;Q]);
else
  error('Wrong selection!');
end

fclose(file);

%if (selection == 2 || selection == 4)
if (false)
  system(sprintf('./%s < %s > %s',programname,infilename,outfilename));
  file = fopen(outfilename,'r');
  T = readTestcase(file);
  fclose(file);
  figure();
  %'Color',[1 1 1],'InvertHardcopy','off','PaperSize',[20.98 29.68]);
  %axes('FontSize',16);
  x = T{1}.m(:,2);
  semilogy(x,T{1}.data(:,2),':','LineWidth',2,'Color',[0,0,0]);
  hold on
  semilogy(x,T{2}.data(:,2),'-.','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{3}.data(:,2),'--','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{4}.data(:,2),'-','LineWidth',2,'Color',[0,0,0]);
  axis([x(1) x(end) 1e-16 1])
  xlabel('N_q');
  ylabel('E_{\infty}','Rotation',0);
end
