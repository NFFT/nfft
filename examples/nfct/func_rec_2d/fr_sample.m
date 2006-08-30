function f = fr_sample( M, min_x, min_y, max_x, max_y, fun, shift_x, shift_y, kind);

if ( kind == ['rand'])
  % random:
  vx = ( max_x - min_x) * rand( M, 1);
  vy = ( max_y - min_y) * rand( M, 1);
  this_M = M;
end

if ( kind == ['grid'])
  % grid:
  M = floor( sqrt( M));
  vy_ = [ min_y + shift_y:( max_y - min_y) / (M-1):max_y + shift_y];
  vy  = [];
  for i = 1:M;
    vy  = [ vy_, vy];
  end;
  vx_ = ones( 1, M)/(M-1);
  vx  = [];
  for i = 1:M;
    vx  = [ vx, shift_x + min_x + ( ( max_x - min_x) * ( i-1) * vx_)];
  end;
  this_M = M * M;
end;

if ( kind == ['data'])
  this_M = M;
  DF = load( fun);
end

fifile = fopen( 'input.dat', 'w');

for i = 1:this_M;

  if ( kind == ['data'])
    
    x = DF( i, 1);
    y = DF( i, 2);
    f = DF( i, 3);
  else
  
    x = vx( i);
    y = vy( i);
    f = feval( fun, x + shift_x, y + shift_y);
  end 
 
  fprintf( fifile, '%f %f %f\n', x, y, f);

end;

fclose( fifile);

