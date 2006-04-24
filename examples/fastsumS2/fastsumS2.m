function T=fastsumS2()

T=[];

% fastsumS2 - Fast summation of radial functions on the sphere
infilename = 'data.in';
outfilenameprefix = 'testcase';
outfilename = 'data.out';
programname = 'fastsumS2';
texfilename = 'table.tex';

selection = menu('fastsumS2 - Fast summation of radial functions on the sphere',...
  'Figure 5.1 (a)','Figure 5.1 (b)','Figure 5.1 (c)','Figure 5.1 (d)','Table 5.2')

file = fopen(infilename,'w');
if (selection > 0 && selection <= 4)
  % Number of testcases
  fprintf(file,'testcases=3\n');
  nodes = [1000, 1000, 1, 0, 1];
  if (selection == 1)
    % Reproduce figure 5.1 (a)
    kernel=0;
    parameters = [0.8];
    m = 4:4:256;
    bound = inline(sprintf('((%d.^(x+1))./(4*pi)).*((2*x+1)./(1-%d)+(2)./((1-%d).^2))',...
      parameters(1),parameters(1),parameters(1)));
  elseif (selection == 2)
    % Reproduce figure 5.1 (b)
    kernel=1;
    parameters = [0.8];
    m = 4:4:256;
    bound = inline(sprintf('((%d.^(x+1))./(4*pi)).*((2*x+1)./(2*(1-%d))+(4*x)./((1-%d).^2)+(4)./((1-%d).^3))',...
      parameters(1),parameters(1),parameters(1),parameters(1)));
  elseif (selection == 3)
    % Reproduce figure 5.1 (c)
    kernel=2;
    parameters = [0.3 7];
    m = 4:4:256;
    bound = inline(sprintf('(1/(pi*sqrt(2*pi))).*(((%d+1).^2)./(%d-0.5)).*(1./((1-%d).^(2*%d+1).*sqrt(sqrt(1-abs(%d))))).*(1./((x-%d).^(%d-0.5)))',...
      parameters(2),parameters(2),parameters(1),parameters(2),parameters(1),parameters(2),parameters(2)));
  elseif (selection == 4)
    % Reproduce figure 5.1 (c)
    kernel=3;
    parameters = [2.5];
    m = 1:1:32;
    bound = inline(sprintf('(sqrt(pi*%d).*(exp(%d)-1).*%d.^(x-0.5))./(gamma(x+0.5))',...
      parameters(1),parameters(1),parameters(1)));
  end
  writeTestcase(file,1,0,0,1,1000,kernel,parameters,m,nodes);
  writeTestcase(file,1,1,3,1,1000,kernel,parameters,m,nodes);
  writeTestcase(file,1,1,6,1,1000,kernel,parameters,m,nodes);
  fclose(file);
  system(sprintf('./%s < %s > %s',programname,infilename,outfilename));
  file = fopen(outfilename,'r');
  T = readTestcase(file);
  fclose(file);
  y4 = feval(bound,T{1}.bandwidths);
  figure('Color',[1 1 1],'InvertHardcopy','off','PaperSize',[20.98 29.68]);
  %axes('FontSize',16);
  semilogy(T{1}.bandwidths,T{1}.data{1}(:,6),'-','LineWidth',2,'Color',[0,0,0]);
  hold on
  x = T{1}.bandwidths;
  semilogy(x,T{2}.data{1}(:,6),'-.','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,T{3}.data{1}(:,6),'--','LineWidth',2,'Color',[0,0,0]);
  semilogy(x,y4,':','LineWidth',2,'Color',[0,0,0]);
  axis([x(1) x(end) 1e-17 1])
  xlabel('M');
  ylabel('E_{\infty}','Rotation',0);
elseif (selection == 5)
  % Number of testcases
  fprintf(file,'testcases=1\n');
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
  kernel=0;
  parameters = [0.6];
  m = 128;
  writeTestcase(file,2,1,6,1,1000,kernel,parameters,m,nodes);
  fclose(file);
  system(sprintf('./%s < %s > %s',programname,infilename,outfilename));
  file = fopen(outfilename,'r');
  T = readTestcase(file);
  fclose(file);
  file = fopen(texfilename,'w');
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
  fclose(file);
else
  error('Invalid selection!');
end

function s = texFormat(d)
if d > 0
  s = sprintf('\\verb#%.1E#',d);
else
  s = '\verb#-#';
end