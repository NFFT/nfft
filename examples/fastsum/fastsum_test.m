% Simple test program for fast NFFT-based summation algorithm.
% Markus Fenn, 2006.

%  d=2;
%  N=4000;
%  M=4000;
%  n=128;
%  m=4;
%  p=3;
%  kernel='multiquadric';
%  c=1/sqrt(N);
%
%  system(sprintf('./fastsum_test %d %d %d %d %d %d %s %e',d,N,M,n,m,p,kernel,c));

N=2000;
M=2000;
kernel='multiquadric';
c=1/sqrt(N);
m=4;
p=3;
n=156;

%random points in circle of radius 0.25-eps_B/2
r=sqrt(rand(N,1))*(0.25-1/32);
phi=rand(N,1)*2*pi;
x=[r.*cos(phi) r.*sin(phi)];

%random coefficients
alpha=rand(N,1)+i*rand(N,1);

%fast NFFT-based summation
f=fastsum(x,alpha,x,kernel,c,m,n,p);