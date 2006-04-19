% Simple test program for fast NFFT-based summation algorithm.
d=2;
N=4000;
M=4000;
n=128;
m=4;
p=3;
kernel='multiquadric';
c=1/sqrt(N);

system(sprintf('./fastsum_test %d %d %d %d %d %d %s %e',d,N,M,n,m,p,kernel,c));
