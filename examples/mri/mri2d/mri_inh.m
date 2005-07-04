N=256;   % points per row / column
M=75520; % 400
%M=75520; % 600
%M=75496; % 900

%load output_phantom_nfft.dat

%M=size(output_phantom_nfft,1)

construct_readout_time( M, 2, 32, 0.00402542373 ); %400
%construct_readout_time( M, 2, 8, 0.00400021193176 );  %900
%construct_readout_time( M, 2, 16,0.00402542372881);%600

%precompute_weights('output_phantom_nfft.dat',M);

disp(2);
system(['./reconstruct_data_inh_3d ' 'output_phantom_nfft.dat ' ...
         int2str(N) ' ' int2str(M)  ' 10 1'])

figure(1);
visualize_data('pics/reallife_400_iter=1',N);
