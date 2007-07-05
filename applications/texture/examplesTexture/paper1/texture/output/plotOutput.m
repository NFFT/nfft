function plotOutput(name, w_hat, option);

newN_values = [0 2.^(1:7)];

close all;
set(0, 'DefaultAxesColorOrder', [0 0 0], ...
		      'DefaultAxesLineStyleOrder', '-+|-.+|--+|:+|-o|-.o|--o|:o|-*');
plot_position = [0.63 0.63 28.45 19.72];
plot_position_save = [3 2 17 15];

files = dir(sprintf('%s.N*.%s', name, w_hat));
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

		data = nan * ones(length(0:2:N), length(newN_values));
	end;	
	
	newN_ind = find(newN == newN_values);

	fid = fopen(curFile, 'r');
	
	results = textscan(fid, '%n%n', 'commentStyle', '#');
	results = results{2};

	if length(results) ~= N+1;
		disp('Inconsistent number of results!');
		disp(length(results));
		disp(N+1);
		return;
	end;	

	data(:, newN_ind) = results(1:2:length(results));

	fclose(fid);
end;	

loglog(1:2:(N+1), data); 

xlabel('l');
ylabel('\epsilon (L, l)');
title(sprintf('L_{max} = %d, N = %d, N'' = %d', N, N1, N2));

for newN_ind = 1:length(newN_values);
	newN = newN_values(newN_ind);
	if newN <= N;
		leg{newN_ind} = sprintf('L = %d', newN);
	end;
end;	

legend(leg, 'location', 'southoutside');

set(gca, 'XTick', newN_values+1);
set(gca, 'XTickLabel', newN_values);
set(gca, 'XMinorTick', 'off');
set(gca, 'YLim', [0.1 * min(min(data)) 10]);
set(gca, 'XLim', [0.9 (N+1)*1.1]);

set(gcf, 'units', 'centimeters');

if strcmp(option, 'print');
	set(gcf, 'paperOrientation', 'landscape');
	set(gcf, 'paperPosition', plot_position);
	print -Plaser1
end;	
if strcmp(option, 'save');
	set(gcf, 'paperOrientation', 'portrait');
	set(gcf, 'paperPosition', plot_position_save);
	saveas(gcf, sprintf('%s.%s.pdf', name, w_hat));
end;	
set(gcf, 'paperOrientation', 'landscape');
set(gcf, 'paperPosition', plot_position);
