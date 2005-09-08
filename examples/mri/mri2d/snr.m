% compute the signal to noise ratio
function [] = snr ( file )
% load original data f
load input_f.dat
% load recontructed data f~
load output_real.dat
load output_imag.dat
output = output_real+i*output_imag;

% compute
result = 20*log10(norm(input_f) / norm(input_f - output));

% write the signal to noise ratio to file
save(file,'result','-ascii');
