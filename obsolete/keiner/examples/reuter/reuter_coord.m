function [theta,phi]=heal_coord(ns)
%
%   Function [theta,phi]=heal_coord(ns);
%
%   This function provide the coordinates of the HealPix pixels.
%
%   Input:            ns       -->    the N_side parameter as explained in the HealPix Primer.
%                                     NB. ns has to be an INTEGER power of
%                                     2.
%                             
%
%   Output            theta    -->    colatitude of the pixels
%                     phi      -->    longitude of the pixels
%              

if log2(ns) ~= fix(log2(ns)), error('NS has to be an INTEGER power of 2'), end

% initialization of the parameters
theta=zeros(12*ns^2,1);
phi=zeros(12*ns^2,1);
ncoo=0;

% North pole
for ii=1:ns-1
    for hh=0:4*ii-1
        ncoo=ncoo+1;
        theta(ncoo)=1-ii^2/(3*ns^2);
        phi(ncoo)=2*pi/(4*ii)*(hh+0.5);
    end
end
ncoo1=ncoo;

% Equator
for ii=ns:3*ns
    for hh=0:4*ns-1
        ncoo=ncoo+1;
        theta(ncoo)=2/(3*ns)*(2*ns-ii);
        phi(ncoo)=2*pi/(4*ns)*(hh+codd(ii));
    end
end

% add the south pole
phi(ncoo+1:end)=phi(ncoo1:-1:1);
theta(ncoo+1:end)=-theta(ncoo1:-1:1);
theta=acos(theta);

return

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function coeff=codd(k);
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if fix(k/2)*2 == k, coeff=0.5;, else, coeff=0;, end

return
