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

N=48;   % points per row / column
Z=48;    % number of slices

% Construct the raw data of the phantom
% and write it to input_f.dat
% To use another example than the phantom
% just put the reshape of your example into input_f.dat
construct_phantom(N,Z);

% Construct the knots in k-space and write them to knots.dat
% Different knots like spiral,rose,radial,radial_3d or linogram can be chosen
% The radial_3d knots just work with construct_data_3d and reconstruct_data_3d
% Then the weights are generated in construct_knots_radial_3d
% M is the number of knots
M = construct_knots_spiral(N,Z);

% First make N^2 1d-FFT, then Z 2d-NFFT on the constructed knots
% and write the output to output_phantom_nfft.dat
if ispc
    cmd='construct_data_2d1d.exe';
else 
    cmd='./construct_data_2d1d';
end
system([cmd ' output_phantom_nfft.dat '...
         int2str(N) ' ' int2str(M) ' ' int2str(Z)]);

% Precompute the weights using voronoi cells
% and write them to weights.dat
precompute_weights_2d('output_phantom_nfft.dat',M,Z);

% First make Z inverse 2d-NFFT, then N^2 inverse 1d-FFT
% and write the output to output_real.dat and output_imag.dat
% The usage is "./reconstruct_data_2d1d filename N M Z ITER WEIGHTS"
% where ITER is the number of iteration and WEIGHTS is 1
% if the weights are used 0 else
if ispc
    cmd='reconstruct_data_2d1d.exe';
else 
    cmd='./reconstruct_data_2d1d';
end
system([cmd ' output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M) ' ' int2str(Z)  ' 3 1']);

% Visualize the three dimensional phantom. Makes a pic of
% every slice and one plot of the N/2 row of the 10th plane.
[~,~] = mkdir('pics_2+1d');
visualize_data('pics_2+1d/pic', N, Z, 1, 'Inverse 2d1d-NFFT - 3. iteration - spiral knots');

% Compute the root mean square
 rms('pics_2+1d/rms.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The same as above but reconstructed with gridding.
% That means first an adjoint 2d-NFFT, then a 1d-FFT.
% The ITER parameter is obsolent and just for compatibility
if ispc
    cmd='reconstruct_data_gridding.exe';
else 
    cmd='./reconstruct_data_gridding';
end
system([cmd ' output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M) ' ' int2str(Z)  ' 0 1']);
[~,~] = mkdir('pics_gridding');
visualize_data('pics_gridding/pic', N, Z, 2, 'Adjoint 2d1d-NFFT (Gridding) - spiral knots');
rms('pics_gridding/rms.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The same as above but reconstructed with a 3d-nfft
if ispc
    cmd='reconstruct_data_3d.exe';
else 
    cmd='./reconstruct_data_3d';
end
system([cmd ' output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M) ' ' int2str(Z)  ' 1 1']);
[~,~] = mkdir('pics_3d');
visualize_data('pics_3d/pic', N, Z, 3, 'Inverse 3d-NFFT - 1. iteration - spiral knots');
rms('pics_3d/rms.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An example with the radial_3d knots
% please uncomment to try

%M = construct_knots_radial_3d(N,N*3/2);

%system(['./construct_data_3d ' 'output_phantom_nfft.dat '...
%         int2str(N) ' ' int2str(M) ' ' int2str(Z)]);

%system(['./reconstruct_data_3d ' 'output_phantom_nfft.dat ' ...
%         int2str(N) ' ' int2str(M) ' ' int2str(Z)  ' 3 1']);
%visualize_data('pics_3d/pic_', N, Z, 4, 'Inverse 3d-NFFT - 3. iteration - radial_3d knots');
%rms('pics_3d/rms.txt');
