N=256;   % points per row / column
M=75520;
arms=16;
%construct_readout_time( M, 2, arms, 0.00402542373 );
%construct_readout_time( M, 2, arms,0.00402542372881);%600


%construct_inh(N);

%out=reshape(phantom(N),1,N*N);
%save input_f.dat -ascii out

%system(['./construct_data_inh_3d ' 'output_phantom_nfft.dat ' ...
%         int2str(N) ' ' int2str(M)])

%precompute_weights('output_phantom_nfft.dat',M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


system(['./reconstruct_data_inh_2d1d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 1 1'])
visualize_data('pics/phantom_2d1d_iter=1',N);
snr('pics/snr_phantom_2d1d_iter=1');

system(['./reconstruct_data_inh_2d1d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 2 1'])
visualize_data('pics/phantom_2d1d_iter=2',N);
snr('pics/snr_phantom_2d1d_iter=2');

system(['./reconstruct_data_inh_2d1d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 5 1'])
visualize_data('pics/phantom_2d1d_iter=5',N);
snr('pics/snr_phantom_2d1d_iter=5');

system(['./reconstruct_data_inh_2d1d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 10 1'])
visualize_data('pics/phantom_2d1d_iter=10',N);
snr('pics/snr_phantom_2d1d_iter=10');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system(['./reconstruct_data_inh_3d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 1 1'])
visualize_data('pics/phantom_3d_iter=1',N);
snr('pics/snr_phantom_3d_iter=1');

system(['./reconstruct_data_inh_3d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 2 1'])
visualize_data('pics/phantom_3d_iter=2',N);
snr('pics/snr_phantom_3d_iter=2');

system(['./reconstruct_data_inh_3d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 5 1'])
visualize_data('pics/phantom_3d_iter=5',N);
snr('pics/snr_phantom_3d_iter=5');

system(['./reconstruct_data_inh_3d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 10 1'])
visualize_data('pics/phantom_3d_iter=10',N);
snr('pics/snr_phantom_3d_iter=10');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system(['./reconstruct_data_2d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 1 1'])
visualize_data('pics/phantom_gridding_iter=1',N);
snr('pics/snr_phantom_gridding_iter=1');

system(['./reconstruct_data_2d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 2 1'])
visualize_data('pics/phantom_gridding_iter=2',N);
snr('pics/snr_phantom_gridding_iter=2');

system(['./reconstruct_data_2d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 5 1'])
visualize_data('pics/phantom_gridding_iter=5',N);
snr('pics/snr_phantom_gridding_iter=5');

system(['./reconstruct_data_2d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 10 1'])
visualize_data('pics/phantom_gridding_iter=10',N);
snr('pics/snr_phantom_gridding_iter=10');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


system(['./reconstruct_data_inh_nnfft ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 1 1'])
visualize_data('pics/phantom_nnfft_iter=1',N);
snr('pics/snr_phantom_nnfft_iter=1');

system(['./reconstruct_data_inh_nnfft ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 2 1'])
visualize_data('pics/phantom_nnfft_iter=2',N);
snr('pics/snr_phantom_nnfft_iter=2');

system(['./reconstruct_data_inh_nnfft ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 5 1'])
visualize_data('pics/phantom_nnfft_iter=5',N);
snr('pics/snr_phantom_nnfft_iter=5');

system(['./reconstruct_data_inh_nnfft ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 10 1'])
visualize_data('pics/phantom_nnfft_iter=10',N);
snr('pics/snr_phantom_nnfft_iter=10');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%