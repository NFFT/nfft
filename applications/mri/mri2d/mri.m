% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%
% This script file is an example of the usage

N=128;   % points per row / column

% Construct the raw data of the phantom
% and write it to input_f.dat
% To use another example than the phantom
% just put the reshape of your example into input_data_f.dat
output = reshape(phantom(N),1,N*N);
save input_f.dat -ascii output

% Construct the knots in k-space and write them to knots.dat
% Different knots like spiral,rose,radial or linogram can be chosen
% M is the number of knots
M = construct_knots_spiral(N,1);

% Make a 2d-NFFT on the constructed knots
% and write the output to output_data_phantom_nfft.dat
if ispc
    cmd='construct_data_2d.exe';
else 
    cmd='./construct_data_2d';
end
system([cmd ' output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)]);

% Precompute the weights using voronoi cells
% and write them to weights.dat
precompute_weights('output_phantom_nfft.dat',M);

% Make an inverse 2d-NFFT
% and write the output to output_real.dat and output_imag.dat
% The usage is "./reconstruct_data_2 filename N M ITER WEIGHTS"
% where ITER is the number of iteration and WEIGHTS is 1
% if the weights are used 0 else
if ispc
    cmd='reconstruct_data_2d.exe';
else 
    cmd='./reconstruct_data_2d';
end
system([cmd ' output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 3 1']);

% Visualize the two dimensional phantom. Make a pic
% and one plot of the N/2 row
[~,~] = mkdir('pics');
visualize_data('pics/pic_2d', N, 1, 'Inverse 2d-NFFT - 3. iteration');

% Compute the root mean square
rms('pics/rms_2d.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The same as above but reconstructed with gridding
% That means an adjoint 2d-NFFT
% The ITER parameter is unused and just for compatibility
if ispc
    cmd='reconstruct_data_gridding.exe';
else 
    cmd='./reconstruct_data_gridding';
end
system([cmd ' output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 5 1']);
visualize_data('pics/pic_gridding', N, 2, '2d-NFFT (Gridding)');
rms('pics/rms_gridding.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

