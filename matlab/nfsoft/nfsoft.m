% This class provides a Matlab interface to the NFSOFT module.
%
% Examples
%   See Matlab scripts simple_test.m.
classdef nfsoft < handle

properties(Dependent=true)
    x;     % nodes (real 3xM matrix)
    fhat;  % Fourier coefficients (complex column vector of length M)
    f;     % samples (complex column vector of length M)
end %properties

properties(SetAccess='private')
    N=[];  % polynomial degree (bandwidth)
    M=[];  % number of nodes
    nfsoft_flags = 0;
    nfft_flags = bitshift(1,12);   % NFFT_OMP_BLOCKWISE_ADJOINT
    nfft_cutoff = 6;
    fpt_kappa = 1000;
    fftw_size;
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

function h=nfsoft(N, M, nfsoft_flags, nfft_flags, nfft_cutoff, fpt_kappa, fftw_size)
% Constructor
%
% plan = nfsoft(N, M)
% plan = nfsoft(N, M, nfsoft_flags)
% plan = nfsoft(N, M, nfsoft_flags, nfft_flags, nfft_cutoff, fpt_kappa, fftw_size)
% 
% INPUT
% N ... polynomial degree (bandwidth)
% M ... number of nodes
% nfsoft_flags  ... can be NFSOFT_NORMALIZED, NFSOFT_REPRESENT,
%   NFSOFT_USE_DPT, NFSOFT_USE_NDFT (default=0)
% nfft_cutoff ... NFFT cutoff parameter (default 6)
% fpt_kappa ... threshold (affects accuracy, 1000 is the default)
% fftw_size ... oversampling (must be >= 2*N+2 and even, default 4*N+4)
% (fftw_size is the size of the fftw transform,
% higher oversampling means better accuracy but more time and memory required)
%
% OUTPUT
%   h   object of class type nfsoft

    h.N=N;
    h.M=M;
    
    if (nargin<2)
        error('Too few arguments')
    end
    if exist('nfsoft_flags','var') && ~isempty(nfsoft_flags)
        h.nfsoft_flags=nfsoft_flags;
    end
    if exist('nfft_flags','var') && ~isempty(nfft_flags)
        h.nfft_flags=nfft_flags;
    end
    if exist('nfft_cutoff','var') && ~isempty(nfft_cutoff)
        h.nfft_cutoff=nfft_cutoff;
    end
    if exist('fpt_kappa','var') && ~isempty(fpt_kappa)
        h.fpt_kappa=fpt_kappa;
    end
    if exist('fftw_size','var') && ~isempty(fftw_size)
        h. fftw_size= fftw_size;
    else
        h. fftw_size=4*h.N+4;
    end
    
    h.plan= nfsoftmex('init',h.N,h.M,h.nfsoft_flags,h.nfft_flags,h.nfft_cutoff,h.fpt_kappa,h.fftw_size);
    h.plan_is_set=true;
end %function

function delete(h)
% Destructor
    if(h.plan_is_set)
         nfsoftmex('finalize',h.plan);
         h.plan_is_set=false;
    end %if
end %function

% Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set.N(h,N)
    if( ~isscalar(N) || (N<1) || round(N)~=N )
        error('The degree N must be a positive integer.')
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

function set.x(h,x)
    if( isempty(x) || ~isnumeric(x) || ~isreal(x) )
        error('The sampling points x have to be real numbers.');
    elseif( size(x,1) ~= 3 || size(x,2) ~= h.M )
        error('The sampling points x must be a 3 x M matrix.');
    elseif( min(x(:))<0 || (max(x(:))>2*pi) || (max(x(2,:))>pi) )
        error('The sampling points x have to be Euler angles on SO(3)');
    else
        nfsoftmex('set_x',h.plan,x);
        nfsoftmex('precompute',h.plan);
        h.x_is_set=true;
    end %if
end %function

function set.fhat(h,fhat)
    nfsoftmex('set_f_hat', h.plan, fhat);
    h.fhat_is_set=true;
end %function

function set.f(h,f)
    if(isempty(f) || ~isnumeric(f) || ~isvector(f) || length(f)~=h.M)
        error('The samples f have to be a vector of length M=%u',h.M);
    else
        nfsoftmex('set_f',h.plan,f(:));
        h.f_is_set=true;
    end %if
end %function

% Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=get.x(h)
    if(h.x_is_set)
        x= nfsoftmex('get_x',h.plan);
    else
        x=[];
    end %if
end %function

function fhat=get.fhat(h)
    if(h.fhat_is_set)
        fhat= nfsoftmex('get_f_hat',h.plan);
    else
        fhat=[];
    end %if
end %funcition

function f=get.f(h)
    if(h.f_is_set)
        f= nfsoftmex('get_f',h.plan);
    else
        f=[];
    end %if
end %function

% User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nfsoft_trafo(h)
% NFSOFT.
%
% nfsoft_trafo(h)
%
% INPUT
%   h  object of class type nfsoft

    if(~h.x_is_set)
        error('Before doing a NFSOFT transform you have to set nodes x.');
    elseif(~h.fhat_is_set)
        error('Before doing a NFSOFT transform you have to set Fourier coefficients in fhat.');
    else
        nfsoftmex('trafo',h.plan);
        h.f_is_set=true;
    end %if
end %function

function nfsoft_adjoint(h)
% Adjoint NFSOFT
%
% nfsoft_adjoint(h)
%
% INPUT
%   h  object of class type nfsoft

    if(~h.x_is_set)
        error('Before doing an adjoint NFSOFT transform you have to set nodes x.');
    elseif(~h.f_is_set)
        error('Before doing an adjoint NFSOFT transform you have to set samples in f.');
    else
        nfsoftmex('adjoint',h.plan);
        h.fhat_is_set=true;
    end %if
end %function

end %methods

end %classdef
