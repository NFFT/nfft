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

% This class provides a Matlab interface to the NFFT module.
%
% Examples
%   See Matlab scripts test_nfft*d.m.
classdef nfft < handle

properties(Dependent=true)
	x;     % nodes (real Mxd matrix)
	fhat;  % fourier coefficients (linearized complex column vector)
	f;     % samples (complex column vector of length M)
	num_threads; % number of threads used for computation
end %properties

properties(SetAccess='private')
	d=[];   % spatial dimension (natural number)
	N=[];   % bandwidth (positive even integer vector)
	M=[];   % number of sampling points (positive integer)
end %properties

properties(Hidden=true,SetAccess='private',GetAccess='private');
	plan=[];                    % nfftmex plan number
	x_is_set=false;             % flag if x is set
	fhat_is_set=false;          % flag if fhat is set
	f_is_set=false;             % flag if f is set
	plan_is_set=false;          % flag if plan was created
	precomputations_done=false; % flag if precomputations were done
end %properties

methods

% Constructer and destructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=nfft(d,N,M,varargin)
% Constructor
%
% h=nfft(1,N,M) for spatial dimension d=1
% h=nfft(2,N,M) for spatial dimension d=2
% h=nfft(3,N,M) for spatial dimension d=3
%
% h=nfft(d,N,M,varargin) for use of nfft_init_guru
% For example
% h=nfft(1,N,M,n,8,bitor(PRE_PHI_HUT,PRE_PSI),FFTW_MEASURE)     for d=1, m=8
% h=nfft(2,N,M,n,n,8,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE) for d=2, m=8
% h=nfft(3,N,M,n,n,n,8,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE) for d=3, m=8
% with n=2^(ceil(log(max(N))/log(2))+1)
% Be careful: There is no error handling with using nfft_init_guru.
% Incorrect inputs can cause a Matlab crash!
%
% INPUT
%   d         spatial dimension
%   N         numbers of nodes in each direction (column vector of length d with positive even integers)
%   M         number of sampling points (positive integer)
%   varargin  parameters for use of nfft_init_guru (see documentation of NFFT for more details)
%
% OUTPUT
%   h   object of class type nfft
	h.d=d;
	if length(N) ~= d
		error('The bandwidth vector N must be an integer vector of length %u for spatial dimension d=%u',d,d);
	end
	h.N=N;
	h.M=M;

	if( 3>nargin )
		error('Wrong number of inputs.');
	elseif( 3==nargin )
		switch d
		case 1
			h.plan=nfftmex('init_1d',N,M);
			h.plan_is_set=true;
		case 2
			h.plan=nfftmex('init_2d',N(1),N(2),M);
			h.plan_is_set=true;
		case 3
			h.plan=nfftmex('init_3d',N(1),N(2),N(3),M);
			h.plan_is_set=true;
		otherwise
			h.plan=nfftmex('init',N,M);
			h.plan_is_set=true;
		end %switch
	else % nfft_init_guru
		args_cell = cell(1,d+2);
		args_cell{1,1} = d;
		for t=1:d
		  args_cell{1,1+t} = N(t);
		end
		args_cell{1,2+d} = M;
		if ~isempty(varargin) && (length(varargin{1}) > 1) % support for n (oversampled N) as vector
			n_array = varargin{1};
			n_array = n_array(:)';
			if length(varargin) > 1
				h.plan=nfftmex('init_guru',[args_cell,num2cell(n_array),varargin{2:end}]);
			else
				h.plan=nfftmex('init_guru',[args_cell,num2cell(n_array)]);
			end
		else % n(1), ..., n(d) (oversampled N) as scalars
			h.plan=nfftmex('init_guru',[args_cell,varargin]);
		end
		h.plan_is_set=true;
	end %if
end %function

function delete(h)
% Destructor
	if(h.plan_is_set)
		nfftmex('finalize',h.plan);
	end %if
end %function

% Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set.d(h,d)
	if( isempty(d) || (d~=round(d)) || (d<1) )
		error('The spatial dimension d has to be a natural number.');
	else
		h.d=d;
	end %if
end %function

function set.N(h,N)
	N = N(:)';
	if( isempty(N) || size(N,1)~=1 || ~isnumeric(N) || ~isreal(N) || (sum(mod(N,2))>0) || ~prod(N>0))
		error('The entries of the bandwidth vector N must be even positive integers.');
	end
	h.N = N;
end

function set.M(h,M)
	if( ~ismatrix(M) || size(M,1)~=1 || size(M,2)~=1)
		error('The number of sampling points M has to be an positive integer.');
	elseif( isempty(M) || ~isnumeric(M) || ~isreal(M) || mod(M,1)~=0 || ~(M>0) )
		error('The number of sampling points M has to be an positive integer.');
	else
		h.M=M;
	end %if
end %function

function set.x(h,x)
	h.x_is_set=false;
	h.precomputations_done=false;

	if( isempty(x) )
		error('The sampling points x have to be real numbers.');
	elseif( ~isnumeric(x) || ~isreal(x) )
		error('The sampling points x have to be real numbers.');
	%elseif( min(x(:))<-1/2 || ~(max(x(:))<1/2) )
	%	error('The sampling points x have to be in the two dimensional Torus [-0.5,0.5)^2');
	elseif( size(x,1)~=h.M || size(x,2)~=h.d )
		error('The sampling points have to be a %ux%u matrix',h.M,h.d);
	end

	x=mod(x+0.5,1)-0.5;
	nfftmex('set_x',h.plan,x.');
	h.x_is_set=true;
	nfftmex('precompute_psi',h.plan)
	h.precomputations_done=true;
end %function

function set.fhat(h,fhat)
	h.fhat_is_set=false;
	n=prod(h.N);

	if( isempty(fhat) || ~isnumeric(fhat))
		error('The Fourier coefficients fhat have to be complex numbers.');
	elseif( size(fhat,1)~=(n) || size(fhat,2)~=1 )
		error('The Fourier coefficients fhat have to be a column vector of length %u.',n);
	end

	if h.d > 1
	  fhat = reshape(fhat, h.N);
	  fhat = permute(fhat, h.d:-1:1);
	  fhat = fhat(:);
	end

	nfftmex('set_f_hat',h.plan,fhat);
	h.fhat_is_set=true;
end %function

function set.f(h,f)
	h.f_is_set=false;
	if(isempty(f) || ~isnumeric(f))
		error('The samples f have to be complex numbers.');
	elseif( size(f,1)~=h.M || size(f,2)~=1 )
		error('The samples f have to be an column vector of length M=%u',h.M);
	else
		nfftmex('set_f',h.plan,f);
		h.f_is_set=true;
	end %if
end %function

function set.num_threads(~,nthreads)
	if nthreads ~= nfft_get_num_threads
		error('Setting the number of threads is not supported');
	end
end %function

% Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=get.x(h)
	if(h.x_is_set)
		x=nfftmex('get_x',h.plan).';
	else
		x=[];
	end %if
end %function

function fhat=get.fhat(h)
	if(h.fhat_is_set)
		fhat=nfftmex('get_f_hat',h.plan);

		if h.d > 1
			fhat = reshape(fhat, h.N(end:-1:1));
			fhat = permute(fhat, h.d:-1:1);
			fhat = fhat(:);
		end
	else
		fhat=[];
	end %if
end %funcition

function f=get.f(h)
	if(h.f_is_set)
		f=nfftmex('get_f',h.plan);
	else
		f=[];
	end %if
end %function

function nthreads=get.num_threads(~)
	nthreads=nfft_get_num_threads;
end %function

% User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nfft_precompute_psi(~)
% Precomputations for NFFT now performed in set.x.
end %function

function ndft_trafo(h)
% NDFT.
%
% ndft_trafo(h)
%
% INPUT
%   h  object of class type nfft

	if(~h.precomputations_done)
		error('Before doing a NFFT transform you have to do precomputations.');
	elseif(~h.fhat_is_set)
		error('Before doing a NFFT transform you have to set Fourier coefficients in fhat.');
	else
		nfftmex('trafo_direct',h.plan);
		h.f_is_set=true;
	end %if
end %function

function nfft_trafo(h)
% NFFT.
%
% nfft_trafo(h)
%
% INPUT
%   h  object of class type nfft

	if(~h.precomputations_done)
		error('Before doing a NFFT transform you have to do precomputations.');
	elseif(~h.fhat_is_set)
		error('Before doing a NFFT transform you have to set Fourier coefficients in fhat.');
	else
		nfftmex('trafo',h.plan);
		h.f_is_set=true;
	end %if
end %function

function ndft_adjoint(h)
% Adjoint NDFT.
%
% ndft_adjoint(h)
%
% INPUT
%   h  object of class type nfft

	if(~h.precomputations_done)
		error('Before doing a adjoint NFFT transform you have to do precomputations.');
	elseif(~h.f_is_set)
		error('Before doing a adjoint NFFT transform you have to set samples in f.');
	else
		nfftmex('adjoint_direct',h.plan);
		h.fhat_is_set=true;
	end %if
end %function

function nfft_adjoint(h)
% Adjoint NFFT
%
% nfft_adjoint(h)
%
% INPUT
%   h  object of class type nfft

	if(~h.precomputations_done)
		error('Before doing a adjoint NFFT transform you have to do precomputations.');
	elseif(~h.f_is_set)
		error('Before doing a adjoint NFFT transform you have to set samples in f.');
	else
		nfftmex('adjoint',h.plan);
		h.fhat_is_set=true;
	end %if
end %function

end %methods

end %classdef

