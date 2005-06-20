% This script file is an example of the usage 

N=128;   % points per row / column

% Construct the raw data of the phantom
% and write it to input_f.dat
% To use another example than the phantom
% just put the reshape of your example into input_data_f.dat
output = reshape(phantom(N),1,N*N);
save input_f.dat -ascii output

% Construct the nodes in k-space and write them to nodes.dat
% Different nodes like spiral,rose,radial or linogram can be chosen
% M is the number of nodes
M = construct_nodes_spiral(N);

% Make a 2d-NFFT on the constructed nodes
% and write the output to output_data_phantom_nfft.dat
system(['./construct_data ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)]);

% Precompute the weights using voronoi cells
% and write them to weights.dat
precompute_weights('output_phantom_nfft.dat',M);

% First make Z inverse 2d-NFFT, then N^2 inverse 1d-FFT
% and write the output to output_real.dat and output_imag.dat
% The usage is "./reconstruct_data_2+1d filename N M ITER WEIGHTS"
% where ITER is the number of iteration and WEIGHTS is 1
% if the weights are used 0 else
system(['./reconstruct_data_2d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 1 1']);

% Visualize the two dimensional phantom. Make a pic
% and one plot of the N/2
visualize_data('pics/pic_2d', N);

% Compute the signal to noise ratio 
%snr('pics/snr_2d.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The same as above but reconstructed with gridding
% That means an adjoint 2d NFFT
% The ITER parameter is obsolent and just for compatibility
%system(['./reconstruct_data_gridding ' 'output_phantom_nfft.dat ' ...
%         int2str(N) ' ' int2str(M)  ' 5 1']);
%visualize_data('pics/pic_gridding', N);
%snr('pics/snr_gridding.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%!rm nodes.dat
%!rm weights.dat
%!rm input_f.dat
%!rm output_phantom_nfft.dat
%!rm output_imag.dat
%!rm output_real.dat

