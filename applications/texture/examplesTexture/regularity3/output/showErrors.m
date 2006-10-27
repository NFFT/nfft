function showErrors(patterns, titles);

%figure();
numplots = size(patterns,1);
lines = floor(sqrt(numplots));
cols = ceil(numplots / lines);
colors = {'k', 'b', 'g', 'r', 'c', 'm', 'y'}
global files;
for j = 1:numplots;
	subplot(lines, cols, j);
    files = {};
    for i = 1:size(patterns,2);
        files{i} = dir(patterns{j,i});
    end;

    for i = 1:length(files);

        c = colors{mod(i-1, length(colors))+1};
		fid = fopen(files{i}(1).name, 'r');

		data = textscan(fid, '%n %n', 'commentStyle', '#');
        nonzero = find(data{2} ~= 0);
        
		semilogy(data{1}(nonzero), data{2}(nonzero), c);
		hold on;
		fclose(fid);
	end;
    title(titles(j));
    legend();
end;
