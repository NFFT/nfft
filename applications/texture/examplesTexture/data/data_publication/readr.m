function r = readr(name);

fid = fopen(name);
textscan(fid, 'Nodes');
textscan(fid, '%*d', 1, 'commentStyle', '#');

r = textscan(fid, '%n %n');
r = [r{1} r{2}];

fclose(fid);