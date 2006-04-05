function fastsumS2()
% fastsumS2 - Fast summation of radial functions on the sphere
infilename = 'data.in';
outfilenameprefix = 'testcase';
outfilename = 'data.out';
programname = 'fastsumS2';

selection = menu('fastsumS2 - Fast summation of radial functions on the sphere',...
  'Figure 5.1 (a)','Figure 5.1 (b)','Figure 5.1 (c)','Figure 5.1 (d)')

file = fopen(infilename,'w');
% Number of testcases
fprintf(file,'testcases=3\n');
nodes = [1000, 1000, 1, 0, 1];
if (selection == 1)
  % Reproduce figure 5.1 (a)
  kernel=0;
  parameters = [0.8];
  m = 4:4:96;
  bound = inline(sprintf('((%d.^(x+1))./(4*pi)).*((2*x+1)./(1-%d)+(2)./((1-%d).^2))',...
    parameters(1),parameters(1),parameters(1)));
end
if (selection == 2)
  % Reproduce figure 5.1 (b)
  kernel=1;
  parameters = [0.8];
  m = 4:4:96;
  bound = inline(sprintf('((%d.^(x+1))./(4*pi)).*((2*x+1)./(2*(1-%d))+(4*x)./((1-%d).^2)+(4)./((1-%d).^3))',...
    parameters(1),parameters(1),parameters(1),parameters(1)));
end
if (selection == 3)
  % Reproduce figure 5.1 (c)
  kernel=3;
  parameters = [0.3 7];
  m = 4:4:128;
  bound = inline(sprintf('(1/(pi*sqrt(2*pi))).*(((%d+1).^2)./(%d-0.5))',...
    parameters(2),parameters(2)));
%  bound = inline(sprintf('(1/(pi*sqrt(2*pi))).*((%d+1).^2./(%d-0.5)).*(1./((1-%d).^(2*%d+1)*sqrt(sqrt(1-abs(%d))))).*(1./((x-%d).^()%d-0.5))',...
%    parameters(2),parameters(2),parameters(1),parameters(2),parameters(1),parameters(2),parameters(2)));
end
if (selection == 4)
  return
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
axes('FontSize',16);
semilogy(T{1}.bandwidths,T{1}.data{1}(:,6),'-','LineWidth',2,'Color',[0,0,0]);
hold on
x = T{1}.bandwidths;
semilogy(x,T{2}.data{1}(:,6),'-.','LineWidth',2,'Color',[0,0,0]);
semilogy(x,T{3}.data{1}(:,6),'--','LineWidth',2,'Color',[0,0,0]);
semilogy(x,y4,':','LineWidth',2,'Color',[0,0,0]);
axis([x(1) x(end) 1e-16 1])
xlabel('M');
ylabel('E_{\infty}','Rotation',0);
