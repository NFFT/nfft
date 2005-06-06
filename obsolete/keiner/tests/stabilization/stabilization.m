f = fopen('stabilization.dat','r');
X = fscanf(f,'%d %d %d %d %f\n',[5,inf])
fclose(f);

thres = unique(X(1,:));

figure;
hold on;
for t = thres
  ind = X(1,:) == t;
  plot(X(2,ind),X(5,ind));
end

xlabel('Bandwidth (M)');
ylabel('Relative number of stabilization steps');
axis tight;
hold off;
print('-depsc',['stabilization.eps']);
