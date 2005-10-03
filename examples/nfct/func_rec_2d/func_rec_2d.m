%    M: evaluate given function at M (random) nodes
%rec_M: reconstruct at rec_M * rec_M grid
function func_rec_2d( M, rec_M)

border_eps=0.0;

% damping weights options
mu=1.2;
c=0.8;

%bandwidth and iterations
N=256;
% for a number of runs use [i0 i1 ... in]
% output is written into seperate files
%iter=[10 20 30 40];
iter=[40];

%misc options (0=false, 1=true)
wait_for_key=0;                  % have matlab wait for input after each run 
show_fig=1;                      % show 3d-figure after computation
make_eps_file=1;                 % create eps-file of figure, show_fig must be set
show_axis=0;                     % have axes in figure

% compile func_rec-2d if neseccary
!make


% create M random 'x y f(x,y)' pairs
% write results to 'input.dat'
% evaluate function 'func.m' at random nodes
% between min_x:max_x and min_y:max_y
% shift x and y by shift_{x,y} (meaning func(x-shift_x, y-shift_y))
% example: fr_sample( M, min_x, max_x, min_y, max_y, 'func', shift_x, shift_y);
% (shift function into first quadrant)

% for now there are "sinx_xysqr: sinc( x^2 + y^2)
%                   "franke:     franke( x ,y)
fr_sample( M, -3.0, -3.0, 3.0, 3.0, 'sinc_xysqr', 3.0, 3.0);


for l = iter

  % start func_rec_2d
  % reads file 'input.dat' and scales input automatically
  % output is written to 'output_l.dat'
  eval( ['!func_rec_2d ' num2str( N) ' ' num2str( M) ' ' num2str( l) ' ' num2str( rec_M) ' ' num2str( border_eps) ' ' num2str( mu) ' ' num2str( c) ' ' ['output_' num2str( l) '.dat'] ]);

  
  % if show_fig is set then view the result
  if( show_fig == 1)

    % load results of reconstruction with l iterations
    A = load( ['output_' num2str( l) '.dat']);

    C = reshape( A(:,3), rec_M, rec_M);
    X = reshape( A(:,1), rec_M, rec_M);
    Y = reshape( A(:,2), rec_M, rec_M);

    fig0 = figure( 1);
    surfl( X, Y, C);

    colormap(bone);                  % set colormap
    shading interp;                  % set shading mode
    
    if( show_axis == 1) axis on; else axis off; end

    view( [ 30, 30]);
    if( make_eps_file == 1) print(fig0, '-deps', ['output_' num2str( l) '.dat.eps']);end  

    if( wait_for_key == 1) pause; end

  end  

end

