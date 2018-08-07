% Source: https://github.com/polarch/Spherical-Harmonic-Transform
% 
% Copyright (c) 2015, Archontis Politis
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% * Neither the name of Spherical-Harmonic-Transform nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [D_l, d_l, R_l] = wignerD(l, alpha, beta, gamma)
%WIGNERD Returns the Wigner D-matrix and d-matrix for SH rotation
%
%   WIGNERD computes the rotation matrices for rotation of functions in the
%   spherical harmonic domain, as given directly by Wigner for complex SHs.
%   Since it much faster to compute the rotation matrices ith recursive
%   formulas, these are included here mostly for comparison.
%
%   Main reference:
%       D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii. 
%       Quantum theory of angular momentum. World Scientific Pub., 1988.
%       p.77 - 4.3(5)
%
%   Inputs: l degree (or band) index
%           alpha, beta, gamma (z, y', z'')-rotation Euler angles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 6/5/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analytic expressions of small-d for order l = 1, 2, included for
% comparison (taken from Varshalovich etal., "Quantum Theory of Angular
% Momentum", 1988, pg.119
% note that the row and column order is inverted from -m:m compared to the 
% book (m:-1:-m)

d1 = @(bta) [   (1+cos(bta))/2      -sin(bta)/sqrt(2)    (1-cos(bta))/2;
                sin(bta)/sqrt(2)    cos(bta)            -sin(bta)/sqrt(2);
                (1-cos(bta))/2      sin(bta)/sqrt(2)   (1+cos(bta))/2];

d2 = @(bta) [   (1+cos(bta))^2/4            -(1+cos(bta))*sin(bta)/2        sqrt(6)*(sin(bta))^2/4          -(1-cos(bta))*sin(bta)/2        (1-cos(bta))^2/4;
                (1+cos(bta))*sin(bta)/2     (1+cos(bta))*(2*cos(bta)-1)/2   -sqrt(3/2)*cos(bta)*sin(bta)    (1-cos(bta))*(2*cos(bta)+1)/2   -(1-cos(bta))*sin(bta)/2;
                sqrt(6)*(sin(bta))^2/4      sqrt(3/2)*sin(bta)*cos(bta)     3*(cos(bta))^2/2-1/2            -sqrt(3/2)*cos(bta)*sin(bta)    sqrt(6)*(sin(bta))^2/4;
                (1-cos(bta))*sin(bta)/2     (1-cos(bta))*(2*cos(bta)+1)/2   sqrt(3/2)*cos(bta)*sin(bta)     (1+cos(bta))*(2*cos(bta)-1)/2   -(1+cos(bta))*sin(bta)/2;
                (1-cos(bta))^2/4            (1-cos(bta))*sin(bta)/2         sqrt(6)*(sin(bta))^2/4          (1+cos(bta))*sin(bta)/2         (1+cos(bta))^2/4];


% precompute factorials
fctrl = factorial(0:2*l).';

% if beta close to zero force small d to identity matrix
D_l = zeros(2*l+1);
if beta <= eps
    d_l = eye(2*l+1);
    D_l = diag(exp(-1i*(alpha+gamma)*(-l:l)));
else
    d_l = zeros(2*l+1);
    for m = -l:l
        for n = -l:l
            % Varshalovich, Eq.4.3.1(5)
            k = (max(0,n-m):min(l+n,l-m)).';
            m1_t = (-1).^k;
            fact_t = sqrt(fctrl(l+n+1)*fctrl(l-n+1)*fctrl(l+m+1)*fctrl(l-m+1)) ./ ...
                        (fctrl(l+n-k+1).*fctrl(l-m-k+1).*fctrl(k+m-n+1).*fctrl(k+1));
            cos_beta = cos(beta/2).^(2*l+n-m-2*k);
            sin_beta = sin(beta/2).^(2*k+m-n);
            d_l_mn = (-1)^(m-n) * sum(m1_t.*fact_t.*cos_beta.*sin_beta);
            d_l(n+l+1,m+l+1) = d_l_mn;
            D_l(n+l+1,m+l+1) = exp(-1i*alpha*m)*d_l_mn*exp(-1i*gamma*n);            
            
        end
    end
end

% Interestingly, if the Wigner-D are used directly they result in an active
% body rotation , instead of rotation of the coordinate system (!). Since
% most operations here are defined in terms of coordinate system rotations,
% the proper result is given by the conjugate of the Wigner-D. Compare for
% example with the rotation matrices using the recursive formula function.
% This is also the same as using e^(ima) * d^l_mn(b) * e^(ing) for the
% entries of the Wigner-D, something that can be found in the literature as
% well as the above formula.

% Compute rotation matrices up to l
R_l = conj(D_l);

end
