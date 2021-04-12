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
function [] = construct_phantom( N, Z )
% constructs a three dimensional Shepp Logan phantom

% The phantom will be a little bit smaller than N
% because of aliasing effects
%N_=N/8*5;
N_=N;

% The specs of the phantom
A = [1.0 0.69 0.92 0.0 0.0 0;
     0.1 0.6624 0.874 0.0 -0.0184 0;
     0.4 0.11 0.31 0.22 0.0 18;
     0.4 0.16 0.41 -0.22 0.0 -18;
     0.55 0.21 0.25 0.0 0.35 0;
     0.55 0.046 0.046 0.0 0.1 0;
     0.55 0.046 0.046 0.0 -0.1 0;
     0.55 0.046 0.023 -0.08 -0.605 0;
     0.55 0.023 0.023 0.0 -0.606 0;
     0.55 0.023 0.046 0.06 -0.605 0];

% compute the phantom
B = zeros(N_,N_,Z);
for y=1:N_
  y_=    y /(N_/2)-1;
  for x=1:N_
    x_=    x /(N_/2)-1;
    r = sqrt(x_*x_+y_*y_);
    if x_==0 && y_==0
      phi=0;
    elseif (x_ > 0) || (x_ >= 0 && y_ >= 0)
        phi=asin(y_/r);
    else
        phi=asin(-y_/r)+pi;
    end
    for z=1:Z
      for l=1:10
          if(((r*cos(pi*A(l,6)/180+phi)+A(l,4))/A(l,2))^2+...
          ((r*sin(pi*A(l,6)/180+phi)+A(l,5))/A(l,3))^2+...
          ((z/(Z/2)-1)/A(l,3))^2 <= 1 )
              B(x,y,z) = B(x,y,z)*(B(x,y,z) > 0.1) .* (l>2) +A(l,1);
          end
      end
    end
  end
end

% place the matrix B in the middle of the matrix C
C=zeros(N,N,Z);
C(N/2-N_/2+1:N/2+N_/2,N/2-N_/2+1:N/2+N_/2,:)=B;
B=C;

% rotate the matrix B
for z=1:Z
  B(:,:,z)=rot90(rot90(rot90((B(:,:,z)))));
end

output=zeros(Z,N*N);

for z_=0:Z-1
  output(z_+1,:)=reshape(B(:,:,z_+1),1,N*N);
end

save input_f.dat -ascii output
