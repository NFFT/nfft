function plotOutput(name, algo, w_hat, measure);

N_values = [5 10 20 40 80];
N1_values = [11 23 41 92 164 308];
N2_values = [128 510 1652 8556 26772 94900];
statusCodes = {'+' '*' 'o'}; 
set(0, 'DefaultAxesColorOrder', [0 0 0], ...
		      'DefaultAxesLineStyleOrder', '-|-.|--|:');
threshold = 10;
testcases = 10;
measure_descr = {'rrn', 'ren', 'rwen'};


data = nan*ones(length(N1_values), length(N_values), testcases);
err = nan*ones(length(N1_values), length(N_values));
status = ones(length(N1_values), length(N_values));


files = dir(sprintf('%s.*.%s.%s.*', name, algo, w_hat));
pattern = sprintf('%s.N_%%d.N1_%%d.N2_%%d.%s.%s.%%d-%%d', name, algo, w_hat);

for count=1:length(files);
	curFile = files(count).name;
	[N, N1, N2, firstcase, lastcase] = strread(curFile, pattern);
	N_ind = find(N == N_values);
	N1_ind = find(N1 == N1_values);
	N2_ind = find(N2 == N2_values);
	if N1_ind ~= N2_ind;
		disp(curFile);
		disp('Inconsistent N1 and N2!');
		disp(N1);
		disp(N2);
		return;
	end;	

	fid = fopen(curFile, 'r');

	results = textscan(fid, '%n%n%n', 'commentStyle', '#');
	results = results{measure};

	if length(results) < lastcase-firstcase+1 + 3;
		disp(curFile);
		disp('Inconsistent number of testcases!');
		disp(length(results)-3);
		disp(lastcase-firstcase+1);
	end;

	data(N1_ind, N_ind, (firstcase+1):(firstcase+length(results)-3)) = ...
																		 results(4:length(results));

	fclose(fid);
end;	
	
for N1_ind = 1:length(N1_values);
	for N_ind = 1:length(N_values);
		ind = find(data(N1_ind, N_ind, :) == data(N1_ind, N_ind, :));
		cur_data = data(N1_ind, N_ind, ind);
		
		if length(cur_data) < testcases;
			status(N1_ind, N_ind) = 3;
		end;	

		if length(cur_data) >= 1;

			err(N1_ind, N_ind) = mean(cur_data);

			if mean(cur_data)/min(cur_data) > threshold || ...
				max(cur_data)/mean(cur_data) > threshold;

				status(N_ind, N1_ind) = 2;
			end;	
		end;
	end;
end;	

figure();

set(gcf, 'paperOrientation', 'landscape');
set(gcf, 'units', 'centimeters');
set(gcf, 'paperPosition', [0.63 0.63 28.45 19.72]);

loglog(N_values', err');

set(gca, 'XTick', [5 10 20 40 80]);
set(gca, 'XTickLabel', ...
	'5 (dim: 286)|10 (dim: 1771)|20 (dim: 12341)|40 (dim: 91881)|80 (dim: 708561)'); 
set(gca, 'XLim', [4 100]);
set(gca, 'YLim', [1e-17 10]);

title(sprintf('Algorithm: %s', algo));
xlabel('L');
ylabel(sprintf('Error Measure: %s', measure_descr{measure}));

for N1_ind = 1:length(N1_values);
	N1 = N1_values(N1_ind);
	N2 = N2_values(N1_ind);
	leg{N1_ind} = sprintf('N = %d, N'' = %d, N*N'' = %d', N1, N2, N1*N2);
end;	
leg_handle = legend(leg, 'location', 'northeastoutside');

set(leg_handle, 'units', 'centimeters');
pos = get(leg_handle, 'position');
lwidth = pos(3);
	
set(gcf, 'units', 'centimeters');	
pos = get(gcf, 'position');
pos(3) = pos(3) + lwidth;
set(gcf, 'position', pos);

for code = 1:3
	ind = find(status ~= code);

	err_code = err;
	err_code(ind) = nan;

	hold on;
	loglog(N_values', err_code', statusCodes{code});
end;
