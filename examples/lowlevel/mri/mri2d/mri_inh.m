% This script file is an example of the usage 

N=128;   % points per row / column
arms=1;  % number of spiral arms

% Construct the spiral knots in k-space and write them to knots.dat
% M is the number of knots
M = construct_knots_spiral(N,arms);

% Construct the file time.dat which contains the readout time
construct_readout_time( M, 2, arms,0.0004);

% Construct a fieldmap and write the output to inh.dat
construct_inh(N);

% Construct the raw data of the phantom
% and write it to input_f.dat
% To use another example than the phantom
% just put the reshape of your example into input_data_f.dat
out=reshape(phantom(N),1,N*N);
save input_f.dat -ascii out

% Construct the k-space data considering the fieldmap
% and write the output to output_data_phantom_nfft.dat
system(['./construct_data_inh_2d1d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)]);

% Precompute the weights using voronoi cells
% and write them to weights.dat
precompute_weights('output_phantom_nfft.dat',M);

% Reconstruct with the 2d1d method
% and write the output to output_real.dat and output_imag.dat
% The usage is "./reconstruct_data_inh_2d1d filename N M ITER WEIGHTS"
% where ITER is the number of iteration and WEIGHTS is 1
% if the weights are used 0 else
% The other methods can be used by replacing 2d1d with 3d or nnfft
system(['./reconstruct_data_inh_2d1d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 3 1']);

% Visualize the two dimensional phantom. Make a pic
% and one plot of the N/2 row
visualize_data('pics/pic_2d1d', N, 1, 'Reconstruction considering the fieldmap');

% Compute the root mean square
rms('pics/rms_2d1d.txt');

% Reconstruct without considering the fieldmap
system(['./reconstruct_data_2d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 3 1']);
visualize_data('pics/pic_2d', N, 2, 'Inverse 2d-NFFT - 3. iteration');
rms('pics/rms_2d.txt');
