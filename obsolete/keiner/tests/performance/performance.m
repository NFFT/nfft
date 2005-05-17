datasets = {'ndsft.dat','Direct NDSFT','ndsft.eps';
            'ndsft_adjoint.dat','Direct adjoint NDSFT','ndsft_adjoint.eps';
            'nfsft.dat','NFSFT','nfsft.eps';
            'nfsft_adjoint.dat','Adjoint NFSFT','nfsft_adjoint.eps'};

for k = 1:size(datasets,1)
  f = fopen(datasets{k,1},'r');
  rows = fscanf(f,'%d\n',1);
  X = zeros(rows,3);
  for n = 1:rows
    X(n,:) = fscanf(f,'%d %d %f\n',3);
  end  
  fclose(f);
  D = unique(X(:,1));
  figure;
  hold on;
  for d=1:length(D)
    index = find(X(:,1)==D(d));
    plot(X(index,2),X(index,3))
  end
  axis tight;
  xlabel('M (bandwidth)');
  ylabel('time/sec');
  title(datasets{k,2});
  hold off;
  print('-depsc',datasets{k,3});
end  
