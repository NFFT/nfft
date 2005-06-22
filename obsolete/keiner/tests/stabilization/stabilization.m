f = fopen('stabilization.dat','r');
X = fscanf(f,'%d %d %d %d %f\n',[5,inf]);
fclose(f);

thres = unique(X(1,:));
maxy = max(X(5,:));
maxx = max(X(2,:));

for t = thres
  figure;
  ind = X(1,:) == t;
  plot(X(2,ind),X(5,ind));
  xlabel('M (bandwidth)');
  ylabel('Relative number of stabilization steps');
  title(['Threshold = ',num2str(t)]);
  axis([0 maxx 0 maxy]);
  print('-depsc',['stabilization',num2str(t),'.eps']);
end
