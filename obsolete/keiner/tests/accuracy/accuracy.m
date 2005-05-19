lsnames = ls('*_*.dat')
pos = find(isspace(lsnames))
filenames = cell(length(pos),1);
filenames{1,1} = lsnames(1:pos(1)-1);
for i = 2:length(pos)
  filenames{i,1} = lsnames(pos(i-1)+1:pos(i)-1);
end
filenames = sort(filenames)
size(filenames);

data = cell(size(filenames,1),4);
thres = [];
ms = [];

for k = 1:size(filenames,1)
  filenames{k}
  f = fopen(filenames{k},'r');
  type = fscanf(f,'%d\n',1);
  mode = fscanf(f,'%d\n',1);
  value = fscanf(f,'%d\n',1);
  rows = fscanf(f,'%d\n',1);
  
  x = zeros(rows,4);
  for n = 1:rows
    x(n,:) = fscanf(f,'%d %f %f %f\n',4);
  end  
  fclose(f);
  
  data{k,1} = type;  
  data{k,2} = mode;  
  data{k,3} = value;  
  data{k,4} = x;
  
  if (mode == 0)
    thres=[thres;value];
  else
    ms=[ms;value];
  end
end  

data

% Plots for different thresholds
thres = unique(thres);
for k = 1:length(thres)
  figure
  hold on;
  for i = 1:size(data,1)
    if (data{i,2} == 0 && data{i,3} == thres(k))
      x = data{i,4};
      if (data{i,1} == 0)
        loglog(x(:,1),x(:,2),'r');
      else
        loglog(x(:,1),x(:,2),'b');      
      end
    end  
  end
  title(['Threshold = ',num2str(thres(k))]);
  xlabel('M (bandwidth)');
  ylabel('Error');
  axis tight;
  hold off;
  print('-depsc',['t',num2str(thres(k)),'eps']);
end


% Plots for different thresholds
ms = unique(ms);
for k = 1:length(ms)
  figure
  hold on;
  for i = 1:size(data,1)
    if (data{i,2} == 1 && data{i,3} == ms(k))
      x = data{i,4};
      if (data{i,1} == 0)
        plot(x(:,1),x(:,2),'r');
      else
        plot(x(:,1),x(:,2),'b');      
      end
    end  
  end
  title(['Bandwidth = ',num2str(ms(k))]);
  xlabel('Threshold');
  ylabel('Error');
  axis tight;
  hold off;
  print('-depsc',['m',num2str(ms(k)),'eps']);
end
