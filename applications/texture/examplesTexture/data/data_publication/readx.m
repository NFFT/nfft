%
% $Id$
function x = readx(name);

fid = fopen(name);
textscan(fid, 'Samples');
textscan(fid, '%*d %*d', 1, 'commentStyle', '#');

x = textscan(fid, '%n + %n');
x = x{1} + x{2};

fclose(fid);
