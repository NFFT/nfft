function plotOutput(name, w_hat, option, max_newN, max_l);

newN_values = [0 2.^(1:7)];

close all;
set(0, 'DefaultAxesColorOrder', [0 0 0], ...
		      'DefaultAxesLineStyleOrder', '-+|--+|:+|-o|--o|:o|-*|--*|:*');
plot_position = [0.63 0.63 28.45 19.72];
plot_position_save = [1 1 15 8];

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

		if nargin < 5;
			max_l = N;
		end;	

		if nargin < 4;
			max_newN = max(newN_values);
		end;

		data = nan * ones(length(0:2:max_l), find(max_newN == newN_values));
	end;	

	if newN > max_newN;
		continue;
	end;	
	
	newN_ind = find(newN == newN_values);

	fid = fopen(curFile, 'r');
	
	results = textscan(fid, '%n%n', 'commentStyle', '#');
	results = results{2};

	if length(results) ~= N+1;
		disp('Inconsistent number of results!');
		disp(length(results));
		disp(N+1);
		continue;
	end;	
	odd_part_deriv = results(2:2:length(results)) - ...
		ones(floor(length(results)/2),1);
	if norm(odd_part_deriv, inf) > 0.1;
		disp('Odd part is not one!')
		newN
		odd_part_deriv
	end;	

	data(:, newN_ind) = results(1:2:(max_l+1));

	fclose(fid);
end;	

ind = find(data < 1e-16)
data(ind) = 1e-16;	

loglog(1:2:(max_l+1), data); 

set(gca, 'fontName', 'Times');
set(gca, 'fontSize', 9);

xlabel('l');
ylabel('\epsilon (L, l)');
title(sprintf('%s, L_0 = %d, N = %d, N'' = %d', name, N, N1, N2));

for newN_ind = 1:length(newN_values);
	newN = newN_values(newN_ind);
	if newN <= max_newN;
		leg{newN_ind} = sprintf('L = %d', newN);
	end;
end;	

%set(gcf, 'resizeFcn', @legResize);

set(gca, 'XTick', newN_values+1);
set(gca, 'XTickLabel', newN_values);
set(gca, 'XMinorTick', 'off');
set(gca, 'XLim', [0.9 (max_l+1)*1.1]);

set(gca, 'YTick', [1e-15 1e-10 1e-5 1]);
set(gca, 'YMinorTick', 'off');
set(gca, 'YLim', [1e-17 10 * max(max(data))]);

set(gcf, 'units', 'centimeters');

legHa = legend(leg, 'location', 'eastoutside');

if nargin > 2;
	if strcmp(option, 'print');
		set(gcf, 'paperOrientation', 'landscape');
		set(gcf, 'paperPosition', plot_position);
		print -Plaser1
	end;	
	if strcmp(option, 'save');
		set(gcf, 'paperOrientation', 'portrait');
		set(gcf, 'paperPosition', plot_position_save);
		saveas(gca, sprintf('%s.pdf', strrep(name, '.', '_')));
	end;	
end;
