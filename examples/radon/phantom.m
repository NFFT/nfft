function I=phantom(N)
% phantom(N)
%  generates the (modified) Shepp-Logan phantom of P. Toft
%  as an NxN matrix.
%
% Reference: Peter Toft: "The Radon Transform - Theory and Implementation", Ph.D. thesis.
%   Department of Mathematical Modelling, Technical University of Denmark, June 1996. 326 pages.

I=zeros(N,N);

k=linspace(-1,1,N); %k=k(end:-1:1);
[x,y]=meshgrid(k);

I = I + ( (x/0.69).^2+(y/0.92).^2 <= 1 ) * 1.0;
I = I + ( (x/0.6624).^2+((y+0.0184)/0.874).^2 <= 1 ) * (-0.8);
I = I + ( (( (cos(-18/360*2*pi)*(x-0.22)+sin(-18/360*2*pi)*y) )/0.11).^2+...
          (( (sin(-18/360*2*pi)*(x-0.22)-cos(-18/360*2*pi)*y) )/0.31).^2 <= 1 ) * (-0.2);
I = I + ( (( (cos(18/360*2*pi)*(x+0.22)+sin(18/360*2*pi)*y) )/0.16).^2+...
          (( (sin(18/360*2*pi)*(x+0.22)-cos(18/360*2*pi)*y) )/0.41).^2 <= 1 ) * (-0.2);
I = I + ( (x/0.21).^2+((y-0.35)/0.25).^2 <= 1 ) * 0.1;
I = I + ( (x/0.046).^2+((y-0.1)/0.046).^2 <= 1 ) * 0.1;
I = I + ( (x/0.046).^2+((y+0.1)/0.046).^2 <= 1 ) * 0.1;
I = I + ( ((x+0.08)/0.046).^2+((y+0.605)/0.023).^2 <= 1 ) * 0.1;
I = I + ( (x/0.023).^2+((y+0.606)/0.023).^2 <= 1 ) * 0.1;
I = I + ( ((x-0.06)/0.023).^2+((y+0.605)/0.046).^2 <= 1 ) * 0.1;

I=flipud(I);
