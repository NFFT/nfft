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
  'Gauss-Legendre','Gauss-Legendre (N=64,128,256, N_Q=4:4:512)','Clenshaw-Curtis','Clenshaw-Curtis (N=64,128,256, N_Q=4:4:512)','HEALPix','Equidistribution')

% Open input data file.
file = fopen(infilename,'w');

if (selection == 1)
  % Set the grid type.
  % 0 = Gauss-Legendre
  %gridtype=0;
  % Set the number of repetitions.
  %repetitions=10;
  % Write the number of testcases.
  %fprintf(file,'testcases=3\n');
  % Write the testcase.
  %writeTestcase(file,1,1,3,1,1000,gridtype,repetitions,[2.^(2:10);2.^(2:10)]);
  %writeTestcase(file,1,1,6,1,1000,gridtype,repetitions,[2.^(2:10);2.^(2:10)]);
  %writeTestcase(file,0,0,0,0,1000,gridtype,repetitions,[2.^(2:6);2.^(2:6)]);
elseif (selection == 2)
  % Set the grid type.
  % 0 = Gauss-Legendre
  gridtype=0;
  % Set the number of repetitions.
  repetitions=1;
  % Set the bandwidhts.
  N=4:4:256;
  % Set the size parameters.
  Q=N;
  % Write the number of testcases.
  fprintf(file,'testcases=4\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,16],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,32],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,64],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,128],repetitions,[N;Q]);
elseif (selection == 3)
  % Set the grid type.
  % 1 = Clenshaw-Curtis
  %gridtype=1;
  % Set the bandwidths
  %m = 4:4:4;
  % Set the number of repetitions. Up to now always 1
  %repetitions=1;
  % Write the number of testcases.
  %fprintf(file,'testcases=1\n');
  % Write the testcase.
  %writeTestcase(file,1,1,3,1,1000,gridtype,repetitions,[m;m]);
  %writeTestcase(file,1,1,6,1,1000,gridtype,repetitions,[m;m]);
  %writeTestcase(file,1,0,6,1,1000,gridtype,repetitions,[m;m]);
elseif (selection == 4)
  % Set the grid type.
  % 1 = Clenshaw-Curtis
  gridtype=1;
  % Set the number of repetitions.
  repetitions=1;
  % Set the bandwidhts.
  N=4:4:256;
  % Set the size parameters.
  Q=N;
  % Write the number of testcases.
  fprintf(file,'testcases=4\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,16],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,32],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,64],repetitions,[N;Q]);
  writeTestcase(file,1,1,6,1,1000,gridtype,[0,128],repetitions,[N;Q]);
elseif (selection == 5)
  % Set the grid type.
  % 2 = HEALPix
  %gridtype=2;
  % Set the bandwidths
  %m = 4:4:8;
  % Set the number of repetitions. Up to now always 1
  %repetitions=5;
  % Write the number of testcases.
  %fprintf(file,'testcases=3\n');
  % Write the testcase.
  %writeTestcase(file,1,1,3,1,1000,gridtype,repetitions,[m;m]);
  %writeTestcase(file,1,1,6,1,1000,gridtype,repetitions,[m;m]);
  %writeTestcase(file,1,0,6,1,1000,gridtype,repetitions,[m;m]);
elseif (selection == 6)
  % Set the grid type.
  % 4 = Equidistribution Example 7.1.11
  %gridtype=3;
  % Set the bandwidths
  %m = 4:4:8;
  % Set the number of repetitions. Up to now always 1
  %repetitions=5;
  % Write the number of testcases.
  %fprintf(file,'testcases=3\n');
  % Write the testcase.
  %writeTestcase(file,1,1,3,1,1000,gridtype,repetitions,m);
  %writeTestcase(file,1,1,6,1,1000,gridtype,repetitions,m);
  %writeTestcase(file,0,1,6,1,1000,gridtype,repetitions,m);
else
  error('Wrong selection!');
end

fclose(file);

if (selection == 2 || selection == 4)
  system(sprintf('./%s < %s > %s',programname,infilename,outfilename));
  file = fopen(outfilename,'r');
  T = readTestcase(file);
  fclose(file);
  figure();
  %'Color',[1 1 1],'InvertHardcopy','off','PaperSize',[20.98 29.68]);
  %axes('FontSize',16);
  x = T{1}.m(:,2);
  semilogy(x,T{1}.data(:,2),'-','LineWidth',2,'Color',[0,0,0]);
  hold on
  semilogy(x,T{2}.data(:,2),'--','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{3}.data(:,2),'-.','LineWidth',2,'Color',[0,0,0]);
  axis([x(1) x(end) 1e-16 1])
  xlabel('N_q');
  ylabel('E_{\infty}','Rotation',0);
end
