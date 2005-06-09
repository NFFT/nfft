% compute the signal to noise ratio
function [] = snr ( file )
% load original data f
load input_f.dat
% load recontructed data f~
load output_real.dat

% compute
result = 10*log10(norm(input_f) / norm(input_f - output_real));

% write the signal to noise ratio to file
save(file,'result','-ascii');
