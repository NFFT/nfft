in=load('output_phantom_nfft2.dat');

out = [in(:,2) in(:,1)];

tmp=in(:,3)+i*in(:,4);
tmp=tmp.*exp(2*pi*i* (-21*in(:,2) + in(:,1)));
out = [out real(tmp) imag(tmp)];

save output_phantom_nfft.dat -ascii out
