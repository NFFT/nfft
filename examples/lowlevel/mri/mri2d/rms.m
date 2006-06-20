% computes the root mean square
function [] = rms ( file )
% load original data f
load input_f.dat
% load recontructed data f~
load output_real.dat
load output_imag.dat
output = output_real+i*output_imag;

% compute
result = (norm(input_f - output) / norm(input_f) );


% write the root mean square to file
save(file,'result','-ascii');
