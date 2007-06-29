function plotAll(name, grid, w_hat);

close all;
algo_values = {'CGNR', 'CGNE'};

for algo_ind = 1:2;
	for measure = 1:3;
		algo = algo_values{algo_ind};
		plotOutput(name, grid, algo, w_hat, measure);
	end;
end;
