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
  'Gauss-Legendre','Clenshaw-Curtis','HEALPix','Equidistribution')

% Open input data file.
file = fopen(infilename,'w');

% Write the number of testcases.
fprintf(file,'testcases=1\n');

% Set the grid type. Up to now always Gauss-Legendre...
% 0 = Gauss-Legendre
% 1 = Clenshaw-Curtis
gridtype=0;

% Set the number of repetitions. Up to now always 1
repetitions=1;

% Set the bandwidths
m = 4:4:96;

% Write the testcase.
writeTestcase(file,1,1,6,1,1000,gridtype,repetitions,m);

fclose(file);

%system(sprintf('./%s < %s > %s',programname,infilename,outfilename));
%file = fopen(outfilename,'r');
%T = readTestcase(file);
%fclose(file);
%y4 = feval(bound,T{1}.bandwidths);
%figure('Color',[1 1 1],'InvertHardcopy','off','PaperSize',[20.98 29.68]);
%axes('FontSize',16);
%semilogy(T{1}.bandwidths,T{1}.data{1}(:,6),'-','LineWidth',2,'Color',[0,0,0]);
%hold on
%x = T{1}.bandwidths;
%semilogy(x,T{2}.data{1}(:,6),'-.','LineWidth',2,'Color',[0,0,0]);
%semilogy(x,T{3}.data{1}(:,6),'--','LineWidth',2,'Color',[0,0,0]);
%semilogy(x,y4,':','LineWidth',2,'Color',[0,0,0]);
%axis([x(1) x(end) 1e-16 1])
%xlabel('M');
%ylabel('E_{\infty}','Rotation',0);
