function showErrors(patterns);

%figure();

for j = 1:size(patterns, 1);
for k = 1:size(patterns, 2);
	subplot(size(patterns, 1), size(patterns, 2), (j-1)*size(patterns, 2) + k);
	files = dir(patterns{j, k})

	for i = 1:length(files);
		files(i).name;
		fid = fopen(files(i).name, 'r');

		for lines = 1:23;
			fgetl(fid);
		end;
	
		data = fscanf(fid, '%d %lg');
		data2 = reshape(data', 2, []);
	
		semilogy(data2(1, :), data2(2, :));
		hold on;
		fclose(fid);
	end;
end;
end;
