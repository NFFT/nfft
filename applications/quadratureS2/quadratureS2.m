function quadratureS2()
%QUADRATURES2
%
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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
% The name of the input file
infilename = 'data.in';
% The name of the output file
outfilename = 'data.out';
% The name of the program
if ispc
    programname='quadratureS2.exe';
else 
    programname='./quadratureS2';
end

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
  S=1024*ones(size(N));
  % Write the number of testcases.
  fprintf(file,'testcases=4\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[3],repetitions,[0],[N;S]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[4],repetitions,[0],[N;S]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[5],repetitions,[0],[N;S]);
  writeTestcase(file,1,1,6,1,1000,0,gridtype,[6],repetitions,[0],[N;S]);
elseif (selection == 3)
  % Set the number of repetitions.
  repetitions=1;
  % Set the size parameters.
  S1 = 16:16:1024;
  S2 = [1,2,4,8,16,32,64,128,256,512];
  %S1 = 16:16:256;
  %S2 = [1,2,4,8,16,32,64,128];
  % Set the bandwidhts.
  N1 = 128*ones(1,length(S1));
  N2 = 128*ones(1,length(S2));
  % Write the number of testcases.
  fprintf(file,'testcases=4\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,0,0,[0,128],repetitions,[0],[N1;S1]);
  writeTestcase(file,1,1,6,1,1000,0,1,[0,128],repetitions,[0],[N1;S1]);
  writeTestcase(file,1,1,6,1,1000,0,2,[0,128],repetitions,[0],[N2;S2]);
  writeTestcase(file,1,1,6,1,1000,0,3,[0,128],repetitions,[0],[N1;S1]);
elseif (selection == 4)
  % Set the grid type.
  % 1 = Clenshaw-Curtis
  gridtype=1;
  % Set the number of repetitions.
  repetitions=1;
  % Set the bandwidhts.
  N=16:16:1024;
  %N=16:16:128;
  % Set the number of nodes.
  Q=N.*N;
  % Set the number of repetitions
  R=ceil(7500./(N.^(1.5))); %horzcat(100:-10:20,20*ones(1,length(N)-length(100:-10:20)));
  % Set the bandwidhts.
  N2=[1,2,4,8,16,32,64,128,256,512,1024];
  %N2=[16,32,64,128];
  % Set the number of nodes.
  Q2=N2.*N2;
  % Set the number of repetitions
  R2=ceil(7500./(N2.^(2))); %horzcat(100:-10:20,20*ones(1,length(N)-length(100:-10:20)));
  % Write the number of testcases.
  fprintf(file,'testcases=2\n');
  % Write the testcases.
  writeTestcase(file,1,1,6,1,1000,1,gridtype,[0,128],repetitions,[0],[N;Q;R]);
  writeTestcase(file,0,1,6,1,1000,1,gridtype,[0,128],repetitions,[0],[N2;Q2;R2]);
else
  error('Wrong selection!');
end

fclose(file);

fprintf('Program in execution. Please be patient! This may take a while!\n');

system(sprintf('%s < %s > %s',programname,infilename,outfilename));
file = fopen(outfilename,'r');
T = readTestcase(file);
fclose(file);
figure('Color',[1 1 1],'InvertHardcopy','off','PaperSize',[20.98 29.68]);
axes('FontSize',16);

if (selection == 1)
%if (false)
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
elseif (selection == 2)
  x = T{1}.parameters(:,1);
  semilogy(x,T{1}.data(:,3),'-.','LineWidth',2,'Color',[0,0,0]);
  hold on
  semilogy(x,T{2}.data(:,3),':','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{3}.data(:,3),'--','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{4}.data(:,3),'-','LineWidth',2,'Color',[0,0,0]);
  axis([x(1) x(end) 1e-12 1e-0])
  xlabel('M');
  ylabel('E_{\infty}','Rotation',0);
elseif (selection == 3)
  x = T{1}.parameters(:,2);
  semilogy(x,T{1}.data(:,3),'-','LineWidth',2,'Color',[0,0,0]);
  hold on
  semilogy(x,T{2}.data(:,3),'--','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{3}.data(:,3),':','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{3}.data(:,3),'.','MarkerSize',2,'Color',[0,0,0]);
  semilogy(x,T{4}.data(:,3),'-.','LineWidth',2,'Color',[0,0,0]);
  axis([x(1) x(end) 1e-12 1e+0])
  xlabel('S');
  ylabel('E_{\infty}','Rotation',0);
elseif (selection == 4)
  x = T{1}.parameters(:,1);
  x2 = T{2}.parameters(:,1);
  semilogy(x,T{1}.data(:,1),'--','LineWidth',2,'Color',[0,0,0]);
  hold on
  semilogy(x2,T{2}.data(:,1),'-','LineWidth',2,'Color',[0,0,0]);
  semilogy(x2,T{2}.data(:,1),'.','MarkerSize',2,'Color',[0,0,0]);
  axis([x(1) x(end) 1e-4 1e+6])
  xlabel('M');
  ylabel('T','Rotation',0);
end
