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

% This class provides a Matlab interface to the NFSFT module.
%
% Examples
%   See Matlab scripts simple_test.m.
classdef nfsft < handle

properties(Dependent=true)
    x;     % nodes (real 2xM matrix)
    fhat;  % Fourier coefficients (complex column vector of length M)
    f;     % samples (complex column vector of length M)
end %properties

properties(SetAccess='private')
    N=[];  % polynomial degree (bandwidth)
    M=[];  % number of nodes
    kappa = 1000;
    nfsft_flags = 0;
    fpt_flags = 0;
    nfft_cutoff = 6;
    nfft_flags = bitshift(1,12);   % NFFT_OMP_BLOCKWISE_ADJOINT
end %properties

properties(Hidden=true,SetAccess='private',GetAccess='private')
    plan=[];
    x_is_set=false;             % flag if x is set
    fhat_is_set=false;          % flag if fhat is set
    f_is_set=false;             % flag if f is set
    plan_is_set=false;          % flag if plan was created
end %properties

methods

% Constructer and destructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=nfsft(N, M, nfsft_flags, kappa, nfft_cutoff, fpt_flags, nfft_flags)
% Constructor
%
% h=nfsft(N,M)
% h=nfsft(N,M,nfsft_flags)
% h=nfsft(N,M,nfsft_flags,kappa,nfft_cutoff,fpt_flags,nfft_flags)
%
% INPUT
%   N           polynomial degree (bandwidth)
%   M           number of sampling points (positive integer)
%   nfsft_flags flags (e.g. NFSFT_NORMALIZED)
%   kappa       threshold (affects accuracy, 1000 is the default)
%   nfft_cutoff NFFT cutoff parameter (default 6)
%
% OUTPUT
%   h   object of class type nfsft

    h.N=N;
    h.M=M;
    
    if (nargin<2)
        error('Too few arguments')
    end
    if exist('nfsft_flags','var') && ~isempty(nfsft_flags)
        h.nfsft_flags=nfsft_flags;
    end
    if exist('kappa','var') && ~isempty(kappa)
        h.kappa=kappa;
    end
    if exist('nfft_cutoff','var') && ~isempty(nfft_cutoff)
        h.nfft_cutoff=nfft_cutoff;
    end
    if exist('fpt_flags','var') && ~isempty(fpt_flags)
        h.fpt_flags=fpt_flags;
    end
    if exist('nfft_flags','var') && ~isempty(nfft_flags)
        h.nfft_flags=nfft_flags;
    end
    
    % Make fhat accessible by the user
    h.nfsft_flags = bitor(h.nfsft_flags,NFSFT_PRESERVE_F_HAT);
    
    nfsftmex('precompute',h.N, h.kappa, h.nfsft_flags, h.fpt_flags);
    h.plan=nfsftmex('init_guru',h.N,h.M,h.nfsft_flags,h.nfft_flags,h.nfft_cutoff);
    h.plan_is_set=true;
    
    % Equispaced nodes are automatically set in nfsft_init
    if bitand(h.nfsft_flags,NFSFT_EQUISPACED)
      h.M = (2*N+2) * (N+2);
      h.x_is_set=true;
    end
end %function

function delete(h)
% Destructor
    if(h.plan_is_set)
        nfsftmex('finalize',h.plan);
        h.plan_is_set=false;
    end %if
end %function

% Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set.N(h,N)
    if( ~isscalar(N) || (N<1) || (N>4096) || round(N)~=N )
        error('The degree N must be a positive integer up to 4096.')
    else
        h.N=N;
    end %if
end %function

function set.M(h,M)
    if( ~isscalar(M) || (M<1) || round(M)~=M )
        error('The number of points M must be a positive integer.')
    else
        h.M=M;
    end %if
end %function

function set.nfsft_flags(h,flags)
    if( ~isscalar(flags) || (flags<0) || round(flags)~=flags )
        error('nfsft_flags must be a nonnegative integer.')
    else
        h.nfsft_flags=flags;
    end %if
end %function

function set.nfft_cutoff(h,c)
    if( ~isscalar(c) || (c<1) || round(c)~=c )
        error('nfft_cutoff must be a positive integer.')
    else
        h.nfft_cutoff=c;
    end %if
end %function

function set.x(h,x)
    if( isempty(x) || ~isnumeric(x) || ~isreal(x) )
        error('The sampling points x have to be real numbers.');
    elseif( size(x,1) ~= 2 || size(x,2) ~= h.M )
        error('The sampling points x must be a 2 x M matrix.');
    elseif( min(x(2,:))<0 || (max(x(2,:))>pi) )
        error('The sampling points x have to be in the two dimensional sphere');
    else
        x(1,:)=mod(x(1,:),2*pi);
        nfsftmex('set_x',h.plan,x);
        h.x_is_set=true;
    end %if
end %function

function set.fhat(h,fhat)
    fhat = f_hat(fhat);
    nfsftmex('set_f_hat_linear', h.plan, fhat.f_hat);
    h.fhat_is_set=true;
end %function

function set.f(h,f)
    if(isempty(f) || ~isnumeric(f))
        error('The samples f have to be complex numbers.');
    elseif( ~isvector(f) || length(f)~=h.M )
        error('The samples f have to be a vector of length M=%u',h.M);
    else
        nfsftmex('set_f',h.plan,f(:));
        h.f_is_set=true;
    end %if
end %function

% Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=get.x(h)
    if(h.x_is_set)
        x=nfsftmex('get_x',h.plan);
    else
        x=[];
    end %if
end %function

function fhat=get.fhat(h)
    if(h.fhat_is_set)
        fhat=f_hat(nfsftmex('get_f_hat_linear',h.plan));
    else
        fhat=f_hat([]);
    end %if
end %funcition

function f=get.f(h)
    if(h.f_is_set)
        f=nfsftmex('get_f',h.plan);
    else
        f=[];
    end %if
end %function

% User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nfsft_trafo(h)
% NFSFT.
%
% nfsft_trafo(h)
%
% INPUT
%   h  object of class type nfft

    if(~h.x_is_set)
        error('Before doing a NFSFT transform you have to set nodes x.');
    elseif(~h.fhat_is_set)
        error('Before doing a NFSFT transform you have to set Fourier coefficients in fhat.');
    else
        nfsftmex('trafo',h.plan);
        h.f_is_set=true;
    end %if
end %function

function nfsft_trafo_direct(h)
% NFSFT.
%
% nfsft_trafo_direct(h)
%
% INPUT
%   h  object of class type nfft

    if(~h.x_is_set)
        error('Before doing a NFSFT transform you have to set nodes x.');
    elseif(~h.fhat_is_set)
        error('Before doing a NFSFT transform you have to set Fourier coefficients in fhat.');
    else
        nfsftmex('trafo_direct',h.plan);
        h.f_is_set=true;
    end %if
end %function

function nfsft_adjoint(h)
% Adjoint NFSFT
%
% nfft_adjoint(h)
%
% INPUT
%   h  object of class type nfsft

    if(~h.x_is_set)
        error('Before doing an adjoint NFSFT transform you have to set nodes x.');
    elseif(~h.f_is_set)
        error('Before doing an adjoint NFSFT transform you have to set samples in f.');
    else
        nfsftmex('adjoint',h.plan);
        h.fhat_is_set=true;
    end %if
end %function

function nfsft_adjoint_direct(h)
% Adjoint NFSFT
%
% nfft_adjoint_direct(h)
%
% INPUT
%   h  object of class type nfsft

    if(~h.x_is_set)
        error('Before doing an adjoint NFSFT transform you have to set nodes x.');
    elseif(~h.f_is_set)
        error('Before doing an adjoint NFSFT transform you have to set samples in f.');
    else
        nfsftmex('adjoint_direct',h.plan);
        h.fhat_is_set=true;
    end %if
end %function

end %methods

end %classdef
