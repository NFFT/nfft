% Copyright (c) 2002, 2020 Jens Keiner, Stefan Kunis, Daniel Potts
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

% This class provides a Matlab interface to the NFCT module.
%
% Examples
%   See Matlab scripts test_nfct*d.m.
classdef nfct < handle

properties(Dependent=true)
	x;     % nodes (real Mxd matrix)
	fhat;  % fourier coefficients (complex column vector of length N1, N1*N2 or N1*N2*N3, for d=2 columnwise linearisation of N1xN2 matrix and for d=3 columnwise linearisation of N1xN2xN3 array)
	f;     % samples (complex column vector of length M)
end %properties

properties(Dependent=true,SetAccess='private');
	N;
end %properties

properties(SetAccess='private')
	d=[];   % spatial dimension (d=1, d=2 or d=3)
	M=[];   % number of sampling points (positive integer)
end %properties

properties(Hidden=true,SetAccess='private',GetAccess='private');
	plan=[];
	N1=[];  % number of nodes in first direction (positive even integer)
	N2=[];  % number of nodes in second direction (positive even integer)
	N3=[];  % number of nodes in third direction (positive even integer)

	x_is_set=false;             % flag if x is set
	fhat_is_set=false;          % flag if fhat is set
	f_is_set=false;             % flag if f is set
	plan_is_set=false;          % flag if plan was created
	precomputations_done=false; % flag if precomputations were done
end %properties

methods

% Constructer and destructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=nfct(d,N,M,varargin)
% Constructor
%
% h=nfct(1,N,M) for spatial dimension d=1
% h=nfct(2,N,M) for spatial dimension d=2
% h=nfct(3,N,M) for spatial dimension d=3
%
% h=nfct(d,N,M,varargin) for use of nfct_init_guru
% For example
% h=nfct(1,N,M,n,7,bitor(PRE_PHI_HUT,PRE_PSI),FFTW_MEASURE)     for d=1, m=7
% h=nfct(2,N,M,n,n,7,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE)   for d=2, m=7
% h=nfct(3,N,M,n,n,n,7,bitor(PRE_PHI_HUT,bitor(PRE_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE) for d=3, m=7
% with n=2^(ceil(log(max(N))/log(2))+1)
% Be careful: There is no error handling with using nfct_init_guru.
% Incorrect inputs can cause a Matlab crash!
%
% INPUT
%   d         spatial dimension (d=1, d=2 or d=3)
%   N         numbers of nodes in each direction (column vector of length d with positive even integers)
%   M         number of sampling points (positive integer)
%   varargin  parameters for use of nfct_init_guru (see documentation of NFCT for more details)
%
% OUTPUT
%   h   object of class type nfct

	h.d=d;

	if( isempty(N) || size(N,1)~=d || size(N,2)~=1)
		error('The numbers of nodes N have to be an integer column vector of length %u for spatial dimension d=%u',d,d);
	else
		h.N=N;
	end %if

	h.M=M;

	if( 3>nargin )
		error('Wrong number of inputs.');
	elseif( 3==nargin )
		switch d
		case 1
			h.plan=nfctmex('init_1d',N,M);
			h.plan_is_set=true;
		case 2
			h.plan=nfctmex('init_2d',N(1),N(2),M);
			h.plan_is_set=true;
		case 3
			h.plan=nfctmex('init_3d',N(1),N(2),N(3),M);
			h.plan_is_set=true;
		otherwise
			error('Invalid spatial dimension d.');
		end %switch
	else % nfct_init_guru
		%disp('You are using nfct_init_guru. This is on your own risk. There will be no error handling. Incorrect inputs can cause a Matlab crash.');
		switch d
		case 1
			args=[{d,N(1),M},varargin];
		case 2
			args=[{d,N(1),N(2),M},varargin];
		case 3
			args=[{d,N(1),N(2),N(3),M},varargin];
		otherwise
			error('Unknown error.');
		end %switch
		h.plan=nfctmex('init_guru',args);
		h.plan_is_set=true;
	end %if
end %function

function delete(h)
% Destructor
	if(h.plan_is_set)
		nfctmex('finalize',h.plan);
		h.plan_is_set=false;
	end %if
end %function

% Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set.d(h,d)
	if( isempty(d) || (d~=1 && d~=2 && d~=3) )
		error('The spatial dimension d has to be d=1, d=2 or d=3.');
	else
		h.d=d;
	end %if
end %function

function set.N(h,N)
	switch length(N)
	case 1
		h.N1=N;
	case 2
		h.N1=N(1);
		h.N2=N(2);
	case 3
		h.N1=N(1);
		h.N2=N(2);
		h.N3=N(3);
	otherwise
		error('Unknown error');
	end %switch
end %function

function set.N1(h,N)
	if( isempty(N) || ~isnumeric(N) || ~isreal(N) || ~(N>0))
		error('The number of the nodes N1 has to be a positive integer.');
	else
		h.N1=N;
	end %if
end %function

function set.N2(h,N)
	if( isempty(N) || ~isnumeric(N) || ~isreal(N) || ~(N>0))
		error('The number of the nodes N2 has to be a positive integer.');
	else
		h.N2=N;
	end %if
end %function

function set.N3(h,N)
	if( isempty(N) || ~isnumeric(N) || ~isreal(N) || ~(N>0))
		error('The number of the nodes N2 has to be a positive integer.');
	else
		h.N3=N;
	end %if
end %function

function set.M(h,M)
	if( ~ismatrix(M) || size(M,1)~=1 || size(M,2)~=1)
		error('The number of sampling pints M has to be an positive integer.');
	elseif( isempty(M) || ~isnumeric(M) || ~isreal(M) || mod(M,1)~=0 || ~(M>0) )
		error('The number of sampling pints M has to be an positive integer.');
	else
		h.M=M;
	end %if
end %function

function set.x(h,x)
	if( isempty(x) )
		error('The sampling points x have to be real numbers.');
	elseif( ~isnumeric(x) || ~isreal(x) )
		error('The sampling points x have to be real numbers.');
	%elseif( min(x(:))<-1/2 || ~(max(x(:))<1/2) )
	%	error('The sampling points x have to be in the two dimensional Torus [-0.5,0.5)^2');
	elseif( size(x,1)~=h.M || size(x,2)~=h.d )
		error('The sampling points have to be a %ux%u matrix',h.M,h.d);
	else
		x=mod(x,0.5);
		nfctmex('set_x',h.plan,x.');
		h.x_is_set=true;
		h.precomputations_done=true; % precompute is called within set_x
	end %if
end %function

function set.fhat(h,fhat)
	switch h.d
	case 1
		n=h.N1;
	case 2
		n=h.N1*h.N2;
	case 3
		n=h.N1*h.N2*h.N3;
	otherwise
		error('Unknown error.');
	end % switch

	if( isempty(fhat) || ~isnumeric(fhat))
		error('The Fourier coefficients fhat have to be complex numbers.');
	elseif( size(fhat,1)~=(n) || size(fhat,2)~=1 )
		error('The Fourier coefficients fhat have to be a column vector of length %u.',n);
	else
		switch h.d
		case 1
			% Do nothing.
		case 2
			% linearization in matlab with column (:) operator is columnwise, in NFCT it is rowwise
			fhat=reshape(fhat,h.N1,h.N2).';
			fhat=fhat(:);
		case 3
			% linearization in matlab with column (:) operator is columnwise, in NFCT it is rowwise
			fhat=reshape(fhat,h.N1,h.N2,h.N3);
			fhat=permute(fhat,[3,2,1]);
			fhat=fhat(:);
		otherwise
			error('Unknown error.');
		end %switch

		nfctmex('set_f_hat',h.plan,fhat);
		h.fhat_is_set=true;
	end %if
end %function

function set.f(h,f)
	if(isempty(f) || ~isnumeric(f))
		error('The samples f have to be complex numbers.');
	elseif( size(f,1)~=h.M || size(f,2)~=1 )
		error('The samples f have to be an column vector of length M=%u',h.M);
	else
		nfctmex('set_f',h.plan,f);
		h.f_is_set=true;
	end %if
end %function

% Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=get.x(h)
	if(h.x_is_set)
		x=nfctmex('get_x',h.plan).';
	else
		x=[];
	end %if
end %function

function fhat=get.fhat(h)
	if(h.fhat_is_set)
		fhat=nfctmex('get_f_hat',h.plan);

		switch h.d
		case 1
			% Do nothing.
		case 2
			% linearization in matlab with column (:) operator is columnwise, in NFCT it is rowwise
			fhat=reshape(fhat,h.N2,h.N1).';
			fhat=fhat(:);
		case 3
			% linearization in matlab with column (:) operator is columnwise, in NFCT it is rowwise
			fhat=reshape(fhat,h.N3,h.N2,h.N1);
			fhat=permute(fhat,[3,2,1]);
			fhat=fhat(:);
		otherwise
			error('Unknown error.');
		end %switch
	else
		fhat=[];
	end %if
end %funcition

function f=get.f(h)
	if(h.f_is_set)
		f=nfctmex('get_f',h.plan);
	else
		f=[];
	end %if
end %function

function N=get.N(h)
	N=[h.N1;h.N2;h.N3];
end %function

% User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function nfct_precompute_psi(h)
% % Precomputations for NFCT.
% 	if(~h.x_is_set)
% 		error('Before doing precomputations you have to set nodes in x.');
% 	else
% 		nfctmex('precompute_psi',h.plan)
% 		h.precomputations_done=true;
% 	end %if
% end %function

function ndct_trafo(h)
% NDCT.
%
% ndct_trafo(h)
%
% INPUT
%   h  object of class type nfct

% 	if(~h.precomputations_done)
% 		error('Before doing a NFCT transform you have to do precomputations.');
	if(~h.x_is_set)
		error('Before doing a NFCT transform you have to set the x nodes.');
	elseif(~h.fhat_is_set)
		error('Before doing a NFCT transform you have to set Fourier coefficients in fhat.');
	else
		nfctmex('trafo_direct',h.plan);
		h.f_is_set=true;
	end %if
end %function

function nfct_trafo(h)
% NFCT.
%
% nfct_trafo(h)
%
% INPUT
%   h  object of class type nfct

% 	if(~h.precomputations_done)
% 		error('Before doing a NFCT transform you have to do precomputations.');
	if(~h.x_is_set)
		error('Before doing a NFCT transform you have to set the x nodes.');
	elseif(~h.fhat_is_set)
		error('Before doing a NFCT transform you have to set Fourier coefficients in fhat.');
	else
		nfctmex('trafo',h.plan);
		h.f_is_set=true;
	end %if
end %function

function ndct_adjoint(h)
% Adjoint NDCT.
%
% ndct_adjoint(h)
%
% INPUT
%   h  object of class type nfct

% 	if(~h.precomputations_done)
% 		error('Before doing a adjoint NFCT transform you have to do precomputations.');
	if(~h.x_is_set)
		error('Before doing an adjoint NFCT transform you have to set the x nodes.');
	elseif(~h.f_is_set)
		error('Before doing an adjoint NFCT transform you have to set samples in f.');
	else
		nfctmex('adjoint_direct',h.plan);
		h.fhat_is_set=true;
	end %if
end %function

function nfct_adjoint(h)
% Adjoint NFCT
%
% nfct_adjoint(h)
%
% INPUT
%   h  object of class type nfct

% 	if(~h.precomputations_done)
% 		error('Before doing a adjoint NFCT transform you have to do precomputations.');
	if(~h.x_is_set)
		error('Before doing an adjoint NFCT transform you have to set the x nodes.');
	elseif(~h.f_is_set)
		error('Before doing a adjoint NFCT transform you have to set samples in f.');
	else
		nfctmex('adjoint',h.plan);
		h.fhat_is_set=true;
	end %if
end %function

end %methods

end %classdef

