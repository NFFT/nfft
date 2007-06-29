function plotOutput(name, w_hat);

newN_values = [0 2.^(0:7)];

close all;
set(0, 'DefaultAxesColorOrder', [0 0 0; 0 0 1], ...
		      'DefaultAxesLineStyleOrder', '-+|-.+|--+|:+|-o|-.o|--o|:o|-*');
plot_position = [0.63 0.63 28.45 19.72];

files = dir(sprintf('%s.*.%s', name, w_hat));
pattern = sprintf('%s.N_%%d.N1_%%d.N2_%%d.newN_%%d.%s', name, w_hat);

for count=1:length(files);
	curFile = files(count).name;
	[curN, curN1, curN2, newN] = strread(curFile, pattern);
	if count > 1;
 		if curN ~= N || curN1 ~= N1 || curN2 ~= N2;

			disp('Inconsistent N, N1, N2!');
		end;	
	else	
		N = curN;
		N1 = curN1;
		N2 = curN2;
	end;	
	
	newN_ind = find(newN == newN_values);

	fid = fopen(curFile, 'r');
	
	results = textscan(fid, '%n%n', 'commentStyle', '#');
	results = results{2};
	results = results(1:2:length(results));

%	hold on;
	loglog(1:2:(N+1), results); 
	leg{count} = sprintf('L = %d', newN);
	
	fclose(fid);
end;	

xlabel('l');
ylabel('\epsilon (L, l)');
title(sprintf('L_{max} = %d, N = %d, N'' = %d', N, N1, N2));
legend(leg, 'location', 'northeastoutside');

set(gca, 'XTick', 1:2:(N+1));
set(gca, 'XTickLabel', 0:2:N);
set(gca, 'XMinorTick', 'off');
set(gca, 'YLim', [1e-17 10]);
set(gca, 'XLim', [0.9 N*1.1]);

set(gcf, 'paperOrientation', 'landscape');
set(gcf, 'units', 'centimeters');
set(gcf, 'paperPosition', plot_position);
