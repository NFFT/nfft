function fastsumS2()
%FASTSUMS2
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
% The input file's name for fastsumS2.c
infilename = 'data.in';
% The output file name.
outfilename = 'data.out';
% The name of the fastsumS2.c executable
if ispc
    programname='fastsumS2.exe';
else 
    programname='./fastsumS2';
end
% The name of the file to write the time measurements table to.
texfilename = 'table.tex';

% Create the menu
selection = menu(...
  ['fastsumS2 - Fast summation of radial functions on the sphere'],...
  'Generate Figure 5.1 (a)','Generate Figure 5.1 (b)',...
  'Generate Figure 5.1 (c)','Generate Figure 5.1 (d)',...
  'Generate Table 5.1')

% Open the file for the input to fastsumS2.c
file = fopen(infilename,'w');

% Check if the error plots are to be generated
if (selection > 0 && selection <= 4)

  % Write the number of testcases to the file.
  fprintf(file,'testcases=3\n');

  % Set the single node set to be used. The values in the vector signify
  % - 1000 Use 1000 source nodes,
  % - 1000 use 1000 target nodes,
  % - 1 Compare with direct evaluation of the sums (1 = Yes, 0 = No),
  % - 0 Don't use precomputed values for direct evaluation (1 = Yes, 0 = No),
  % - 1 Do the summation process 1 time.
  nodes = [1000, 1000, 1, 0, 1];

  % Switch by kernel type.
  if (selection == 1)

    % Generate figure 5.1 (a)

    % Set the kernel to be the Abel-Poisson kernel.
    kernel = 0;
    % Set the parameter h.
    parameters = [0.8];
    % Set the cut-off bandwidth
    m = 4:4:256;

    % Generate the theoretical error bound function from [1].
    bound = inline(sprintf(['((%d.^(x+1))./(4*pi)).*((2*x+1)./(1-%d)+(2)./',...
        '((1-%d).^2))'],parameters(1),parameters(1),parameters(1)));

  elseif (selection == 2)

    % Generate figure 5.1 (b)

    % Set the kernel to be the singularity kernel.
    kernel=1;

    % Set the parameter h.
    parameters = [0.8];

    % Set the cut-off bandwidth
    m = 4:4:256;

    % Generate the theoretical error bound function from [1].
    bound = inline(sprintf(['((%d.^(x+1))./(4*pi)).*((2*x+1)./(2*(1-%d))+',...
      '(4*x)./((1-%d).^2)+(4)./((1-%d).^3))'], parameters(1),...
      parameters(1),parameters(1),parameters(1)));

  elseif (selection == 3)

    % Generate figure 5.1 (c)

    % Set the kernel to be the locally supported kernel.
    kernel=2;

    % Set the parameter h and lambda.
    parameters = [0.3 7];

    % Set the cut-off bandwidth
    m = 4:4:256;

    % Generate the theoretical error bound function from [1].
    bound = inline(sprintf(['(1/(pi*sqrt(2*pi))).*(((%d+1).^2)./(%d-0.5)).*',...
      '(1./((1-%d).^(2*%d+1).*sqrt(sqrt(1-abs(%d))))).*(1./((x-%d).^',...
      '(%d-0.5)))'], parameters(2), parameters(2), parameters(1),...
      parameters(2), parameters(1), parameters(2), parameters(2)));

  elseif (selection == 4)

    % Generate figure 5.1 (c)

    % Set the kernel to be the spherical gaussian kernel.
    kernel=3;

    % Set the parameter sigma.
    parameters = [2.5];

    % Set the cut-off bandwidth
    m = 1:1:32;

    % Generate the theoretical error bound function from [1].
    bound = inline(sprintf('(sqrt(pi*%d).*(exp(%d)-1).*%d.^(x-0.5))./(gamma(x+0.5))',...
      parameters(1),parameters(1),parameters(1)));

  end

  % Write three testcases for testing three different algorithms and parameters
  % to the input file.

  % NFSFT algorithm with direct NDFT and direct DPT
  writeTestcase(file,1,0,0,1,1000,kernel,parameters,m,nodes);

  % NFSFT algorithm with NFFT (cut-off parameter 3) and FPT algorithm
  % (threshold 1000.0)
  writeTestcase(file,1,1,3,1,1000,kernel,parameters,m,nodes);

  % NFSFT algorithm with NFFT (cut-off parameter 6) and FPT algorithm
  % (threshold 1000.0)
  writeTestcase(file,1,1,6,1,1000,kernel,parameters,m,nodes);

  % Close the input file.
  fclose(file);

  % Call fastsumS2.c with the generated input file writing the output to the
  % output file.
  system(sprintf('%s < %s > %s',programname,infilename,outfilename));

  % Open the output file.
  file = fopen(outfilename,'r');

  % Read the testcases into a cell array.
  T = readTestcase(file);

  % Close the output file.
  fclose(file);

  % Generate the values of the theoretical error bound function.
  y4 = feval(bound,T{1}.bandwidths);

  % Generate a new figure.
  figure('Color',[1 1 1],'InvertHardcopy','off','PaperSize',[20.98 29.68]);

  % Add the different error curves to the figure

  % Get the cut-off degree values.
  x = T{1}.bandwidths;

  % First testcase
  semilogy(x,T{1}.data{1}(:,6),'-','LineWidth',2,'Color',[0,0,0]);

  % Prevent old curves from being deleted by adding the next curve
  hold on

  % Second testcase
  semilogy(x,T{2}.data{1}(:,6),'-.','LineWidth',2,'Color',[0,0,0]);

  % Third testcase
  semilogy(x,T{3}.data{1}(:,6),'--','LineWidth',2,'Color',[0,0,0]);

  % Theoretical error bound
  semilogy(x,y4,':','LineWidth',2,'Color',[0,0,0]);

  % Adjust the axis limits.
  axis([x(1) x(end) 1e-17 1])

  % Add axis labels.
  xlabel('M');
  ylabel('E_{\infty}','Rotation',0);

elseif (selection == 5)

  % Generate Table 5.1

  % Write the number of testcases to the file.
  fprintf(file,'testcases=1\n');

  % Set node sets to be used. We use L = D source and target nodes as a power
  % of two from 2^6 up to 2^21. Up to L = D = 2^18, we use the direct sum
  % evaluation and compute the error E_infty. Up to L = D = 2^12 we also use
  % precomputed kernel function values to compare the time needed. For small
  % node numbers and computation times, we use repetitions to obtain time
  % measurements averaged over multiple computations.
  nodes = [...
     2^6,  2^6, 1, 1, 60;...
     2^7,  2^7, 1, 1, 50;...
     2^8,  2^8, 1, 1, 35;...
     2^9,  2^9, 1, 1, 20;...
    2^10, 2^10, 1, 1, 10;...
    2^11, 2^11, 1, 1, 10;...
    2^12, 2^12, 1, 1, 10;...
    2^13, 2^13, 1, 0, 1;...
    2^14, 2^14, 1, 0, 1;...
    2^15, 2^15, 1, 0, 1;...
    2^16, 2^16, 1, 0, 1;...
    2^17, 2^17, 1, 0, 1;...
    2^18, 2^18, 1, 0, 1;...
    2^19, 2^19, 0, 0, 1;...
    2^20, 2^20, 0, 0, 1;...
    2^21, 2^21, 0, 0, 1];

  % Set the kernel to be the Abel-Poisson kernel.
  kernel = 0;

  % Set the parameter h.
  parameters = [0.6];

  %Set the cut-off degree.
  m = 128;

  % Write the testcase to the input file.
  writeTestcase(file,2,1,6,1,1000,kernel,parameters,m,nodes);

  % Close the input file.
  fclose(file);

  % Call fastsumS2.c with the generated input file writing the output to the
  % output file.
  system(sprintf('%s < %s > %s',programname,infilename,outfilename));

  % Open the output file.
  file = fopen(outfilename,'r');

  % Read the testcases into a cell array.
  T = readTestcase(file);

  % Close the output file.
  fclose(file);

  % Open the file to write the table in TeX format to.
  file = fopen(texfilename,'w');

  % Generate the table in TeX format and write it to the output file.
  fprintf(file,'\\begin{table}[ht!]\n  \\begin{center}\n    ');
  fprintf(file,'\\begin{tabular}{r|r|r|r|r|r}\n      ');
  fprintf(file,'$L = D$ &      direct alg.   &     w/pre-comp.    &      ');
  fprintf(file,'FS, NDSFT &      FS, NFSFT & error ');
  fprintf(file,'$E_{\\infty}$\\\\\\hline\\\\[-2.0ex]\n');
  for i = 1:size(T{1}.data,2)
    fprintf(file,'      $2^{%d}$ & %s & %s & %s & %s & %s\\\\\n',...
      round(log2(T{1}.nodes(i,1))),texFormat(T{1}.data{i}(1)),...
      texFormat(T{1}.data{i}(2)),texFormat(T{1}.data{i}(3)),...
      texFormat(T{1}.data{i}(4)),texFormat(T{1}.data{i}(6)));
  end
  fprintf(file,'    \\end{tabular}\n');
  fprintf(file,'  \\end{center}\n');
  fprintf(file,'\\end{table}\n');

  % Close the output file.
  fclose(file);
else
  % Error due to invalid selection.
  error('Invalid selection!');
end

% End of the function
return;


function s = texFormat(d)
% TEXFORMAT Output positive numbers in TeX verbatim font style
%   TEXFORMAT(d) converts any number d > 0 to a string in TeX verbatim
%   font style. The number is converted into scientfic format with one decimal
%   digit precision. If d <= 0 a dash in Tex verbatim font style is returned in
%   place of the number.

% Check if number is greater zero.
if d > 0

  % Convert number to TeX verbatim font style and scientif format.
  s = sprintf('\\verb#%.1E#',d);

else

  % Return a dash in TeX verbatim font style.
  s = '\verb#-#';

end
