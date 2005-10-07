function f = fr_sample( M, min_x, min_y, max_x, max_y, fun, shift_x, shift_y);

vx = ( max_x - min_x) * rand( M, 1);
vy = ( max_y - min_y) * rand( M, 1);

ffile = fopen( 'input.dat', 'w');

for i = 1:M;

  x = vx( i) - shift_x;
  y = vy( i) - shift_y;
  
  f = feval( fun, x, y);

  fprintf( ffile, '%f %f %f\n', vx( i), vy( i), f);

end;

fclose(ffile);

