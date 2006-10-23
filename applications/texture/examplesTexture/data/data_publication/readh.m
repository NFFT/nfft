function h = readh(name);

fid = fopen(name);
textscan(fid, 'Polefigures');
textscan(fid, '%*d', 1, 'commentStyle', '#');

h = textscan(fid, '%n %n');
h = [h{1} h{2}];

fclose(fid);