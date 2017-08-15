function plotGrid(gridType,p)
%plotGrid - Plot spherical quadrature grids
%
%   plotGrid(gridtype) plots a quadrature grid. The parameter gridType specifies
%   the point set to plot and can take the following values:
%   0: Gauss-Legendre quadrature grid,
%   1: Clenshaw-Curtis qudrature grid,
%   2: HEALPix point set,
%   3: Equidistribution point set.
%
%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts

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
%
if (gridType == 0)
  [theta,w] = lgwt(p(1)+1,-1,1);
  theta = acos(theta);
  phi = 2*pi*(0:(2*p(1)+1))./(2*p(1)+2);
  [Y,X] = meshgrid(theta,phi);
  gridName = 'Gauss-Legendre';
  temp = repmat(w,[1,length(phi)])';
  [Y(:),X(:),temp(:)]
  size(Y(:))
elseif (gridType == 1)
  theta = pi*(0:2*p(1))./(2*p(1));
  phi = 2*pi*(0:(2*p(1)+1))./(2*p(1)+2);
  [X,Y] = meshgrid(phi,theta);
  gridName = 'Clenshaw-Curtis';

elseif (gridType == 2)
  if log2(p(1)) ~= fix(log2(p(1)))
    error('NS has to be an INTEGER power of 2')
  end

  % initialization of the parameters
  Y=zeros(12*p(1)^2,1);
  X=zeros(12*p(1)^2,1);
  ncoo=0;

  % North pole
  for ii=1:p(1)-1
    for hh=0:4*ii-1
      ncoo=ncoo+1;
      Y(ncoo)=1-ii^2/(3*p(1)^2);
      X(ncoo)=2*pi/(4*ii)*(hh+0.5);
    end
  end
  ncoo1=ncoo;

  % Equator
  for ii=p(1):3*p(1)
    for hh=0:4*p(1)-1
      ncoo=ncoo+1;
      Y(ncoo)=2/(3*p(1))*(2*p(1)-ii);
      X(ncoo)=2*pi/(4*p(1))*(hh+codd(ii));
    end
  end

  % add the south pole
  X(ncoo+1:end)=X(ncoo1:-1:1);
  Y(ncoo+1:end)=-Y(ncoo1:-1:1);
  Y=acos(Y);

  gridName = 'HEALPix';

elseif (gridType == 3)

  X = zeros(2+4*(floor((p(1)+1)/2))*floor(p(1)/2),1);
  Y = zeros(2+4*(floor((p(1)+1)/2))*floor(p(1)/2),1);
  X(1) = 0.0;
  Y(1) = 0.0;
  d = 2;
  for k = 1:p(1)-1
    thetai = (k*pi)/p(1);
    if (k<=(p(1)/2))
      gammai = 4*k;
    else
      gammai = 4*(p(1)-k);
    end
    for n = 1:gammai
      Y(d) = thetai;
      X(d) = ((n-0.5)*((2.0*pi)/gammai));
      d = d+1;
    end
  end
  Y(d) = pi;
  X(d) = 0.0;
  size(X)

  gridName = 'Equidistribution Example 7.1.11';

else
  error('Wrong grid type!');
end

figure;
plot(X,Y,'ko');
axis equal;
axis([-0.1 2*pi+0.1 -0.1 pi+0.1]);
xlabel('Phi');
ylabel('Theta');
title(gridName);
return;

function coeff=codd(k);
  if fix(k/2)*2 == k, coeff=0.5;, else, coeff=0;, end
return
