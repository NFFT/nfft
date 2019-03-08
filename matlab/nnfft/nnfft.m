
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

% This class provides a Matlab interface to the NNFFT module.
%
% Examples
%   See Matlab scripts test_nnfft*d.m.
classdef nnfft < handle

properties(Dependent=true)
	x;     % nodes in time/spatial domain (real Mxd matrix)
	v;		% nodes in fourier domain (real N_totalxd matrix)
	fhat;  % fourier coefficients (complex column vector of length N_total)
	f;     % samples (complex column vector of length M)
end %properties

properties(Dependent=true,SetAccess='private');
	N; %cut-off-frequencies
	N1; %sigma*N
end %properties

properties(SetAccess='private')
	d=[];   % spatial dimension (d=1, d=2 or d=3)
	N_total=[];
	M=[];   % number of sampling points (positive integer)
end %properties

properties(Hidden=true,SetAccess='private',GetAccess='private');
	plan=[];
	N_1=[];  % number of nodes in first direction (positive even integer)
	N_2=[];  % number of nodes in second direction (positive even integer)
	N_3=[];  % number of nodes in third direction (positive even integer)
	N1_1=[]; % number of nodes in first direction (positive even integer)
	N1_2=[]; % number of nodes in second direction (positive even integer)
	N1_3=[]; % number of nodes in third direction (positive even integer)

	x_is_set=false;             % flag if x is set
	v_is_set=false;			    % flag if v is set
	fhat_is_set=false;          % flag if fhat is set
	f_is_set=false;             % flag if f is set
	plan_is_set=false;          % flag if plan was created
	precomputations_done=false; % flag if precomputations were done
end %properties

methods

% Constructer and destructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=nnfft(d,N_total,M,N,varargin)
% Constructor
%
% h=nnfft(1,N_total,M,N) for spatial dimension d=1
% h=nnfft(2,N_total,M,N) for spatial dimension d=2
% h=nnfft(3,N_total,M,N) for spatial dimension d=3
%
% h=nnfft(d,N_total,M,N,varargin) for use of nnfft_init_guru
% For example
% h=nnfft(1,N_total,M,N,N1,7,bitor(PRE_PHI_HUT,PRE_PSI))     for d=1, m=7
% h=nnfft(2,N_total,M,N,N1_1,N1_2,7,bitor(PRE_PHI_HUT,PRE_PSI)) for d=2, m=7
% h=nnfft(3,N_total,M,N,N1_1,N1_2,N1_3,7,bitor(PRE_PHI_HUT,PRE_PSI)) for d=3, m=7
% with N1=sigma*N   ; n=N1
% NOT IMPLEMENTED: Be careful: There is no error handling with using nfft_init_guru.
% Incorrect inputs can cause a Matlab crash!
%
% INPUT
%   d         spatial dimension (d=1, d=2 or d=3)
%	N_total  Total number of Fourier coefficients.
%   M         number of sampling points (positive integer)
%   N			cutt-off-frequencies, numbers of nodes in each direction (column vector of length d with positive even integers)
%   varargin  parameters for use of nnfft_init_guru (see documentation of NFFT for more details)
%
% OUTPUT
%   h   object of class type nnfft

	h.d=d;
	h.N_total=N_total;
	h.M=M;
	if( isempty(N) || size(N,1)~=d || size(N,2)~=1)
		error('The numbers of nodes N have to be an integer column vector of length %u for spatial dimension d=%u',d,d);
	else
		h.N=N;
	end %if

	if( 4>nargin )
		error('Wrong number of inputs.');
	elseif( 4==nargin )
		
		switch d
		case 1
			h.plan=nnfftmex('init_1d',N_total,M,N);
			h.plan_is_set=true;
		case 2
			h.plan=nnfftmex('init_2d',N_total,M,N(1),N(2));
			h.plan_is_set=true;
		case 3
			h.plan=nnfftmex('init_3d',N_total,M,N(1),N(2),N(3));
			h.plan_is_set=true;
		otherwise
			error('Invalid spatial dimension d.');
		end %switch
	else % nnfft_init_guru
		disp('You are using nnfft_init_guru. This is on your own risk. There will be no error handling. Incorrect inputs can cause a Matlab crash.');
		
		switch d
		case 1
			args=[{d,N_total,M,N(1)},varargin];
		case 2
			%disp(N_total);disp(N);
			args=[{d,N_total,M,N(1),N(2)},varargin];
		case 3
			args=[{d,N_total,M,N(1),N(2),N(3)},varargin];
		otherwise
			error('Unknown error.');
		end %switch
		h.plan=nnfftmex('init_guru',args);
		h.plan_is_set=true;
	end %if
end %function

function delete(h)
% Destructor
	if(h.plan_is_set)
		nnfftmex('finalize',h.plan);
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

function set.N_total(h,k)
	if( ndims(k)~=2 || size(k,1)~=1 || size(k,2)~=1)
		error('The number of sampling pints N_total has to be an positive integer.');
	elseif( isempty(k) || ~isnumeric(k) || ~isreal(k) || mod(k,1)~=0 || ~(k>0) )
		error('The number of sampling pints N_total has to be an positive integer.');
	else
		h.N_total=k;
	end %if
end %function

function set.N(h,N)
	switch length(N)
	case 1
		h.N_1=N;
	case 2
		h.N_1=N(1);
		h.N_2=N(2);
	case 3
		h.N_1=N(1);
		h.N_2=N(2);
		h.N_3=N(3);
	otherwise
		error('Unknown error');
	end %switch
end %function

function set.N_1(h,N_1)
	if( isempty(N_1) || ~isnumeric(N_1) || ~isreal(N_1) || (mod(N_1,2)~=0) || ~(N_1>0))
		error('The number of the nodes N_1 has to be an even positive integer.');
	else
		h.N_1=N_1;
	end %if
end %function

function set.N_2(h,N_2)
	if( isempty(N_2) || ~isnumeric(N_2) || ~isreal(N_2) || (mod(N_2,2)~=0) || ~(N_2>0))
		error('The number of the nodes N_2 has to be an even positive integer.');
	else
		h.N_2=N_2;
	end %if
end %function

function set.N_3(h,N_3)
	if( isempty(N_3) || ~isnumeric(N_3) || ~isreal(N_3) || (mod(N_3,2)~=0) || ~(N_3>0))
		error('The number of the nodes N_3 has to be an even positive integer.');
	else
		h.N_3=N_3;
	end %if
end %function


function set.N1(h,N1)
	switch length(N1)
	case 1
		h.N1_1=N1;
	case 2
		h.N1_1=N1(1);
		h.N1_2=N1(2);
	case 3
		h.N1_1=N1(1);
		h.N1_2=N1(2);
		h.N1_3=N1(3);
	otherwise
		error('Unknown error');
	end %switch
end %function

function set.N1_1(h,N1_1)
	if( isempty(N1_1) || ~isnumeric(N1_1) || ~isreal(N1_1) || (mod(N1_1,2)~=0) || ~(N1_1>0))
		error('The number of the nodes N1_1 has to be an even positive integer.');
	else
		h.N1_1=N1_1;
	end %if
end %function

function set.N1_2(h,N1_2)
	if( isempty(N1_2) || ~isnumeric(N1_2) || ~isreal(N1_2) || (mod(N1_2,2)~=0) || ~(N1_2>0))
		error('The number of the nodes N1_2 has to be an even positive integer.');
	else
		h.N1_2=N1_2;
	end %if
end %function

function set.N1_3(h,N1_3)
	if( isempty(N1_3) || ~isnumeric(N1_3) || ~isreal(N1_3) || (mod(N1_3,2)~=0) || ~(N1_3>0))
		error('The number of the nodes N1_3 has to be an even positive integer.');
	else
		h.N1_3=N1_3;
	end %if
end %function

function set.M(h,M)
	if( ndims(M)~=2 || size(M,1)~=1 || size(M,2)~=1)
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
		x=mod(x+0.5,1)-0.5;
		nnfftmex('set_x',h.plan,x.');
		h.x_is_set=true;
		h.precomputations_done=false;
	end %if
end %function

function set.v(h,v)
	if( isempty(v) )
		error('The sampling points v for fourier domain have to be real numbers.');
	elseif( ~isnumeric(v) || ~isreal(v) )
		error('The sampling points v for fourier domain have to be real numbers.');
	%elseif( min(v(:))<-1/2 || ~(max(v(:))<1/2) )
	%	error('The sampling points v for fourier domain have to be in the two dimensional Torus [-0.5,0.5)^2');
	elseif( size(v,1)~=h.N_total || size(v,2)~=h.d )
		error('The sampling points for fourier domain have to be a %uv%u matrix',h.N_total,h.d);
	else
		v=mod(v+0.5,1)-0.5;
		nnfftmex('set_v',h.plan,v.');
		h.v_is_set=true;
		h.precomputations_done=false;
	end %if
end %function

function set.fhat(h,fhat)
	%switch h.d
	%case 1
		%n=h.N_1;
	%case 2
	%	n=h.N_1*h.N_2;
	%case 3
	%	n=h.N_1*h.N_2*h.N_3;
	%otherwise
	%	error('Unknown error.');
	%end % switch

	n=h.N_total;

	if( isempty(fhat) || ~isnumeric(fhat))
		error('The Fourier coefficients fhat have to be complex numbers.');
	elseif( size(fhat,1)~=(n) || size(fhat,2)~=1 )
		error('The Fourier coefficients fhat have to be a column vector of length %u.',n);
	else
		nnfftmex('set_f_hat',h.plan,fhat);
		h.fhat_is_set=true;
	end %if
end %function

function set.f(h,f)
	if(isempty(f) || ~isnumeric(f))
		error('The samples f have to be complex numbers.');
	elseif( size(f,1)~=h.M || size(f,2)~=1 )
		error('The samples f have to be an column vector of length M=%u',h.M);
	else
		nnfftmex('set_f',h.plan,f);
		h.f_is_set=true;
	end %if
end %function

% Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=get.x(h)
	if(h.x_is_set)
		x=nnfftmex('get_x',h.plan).';
	else
		x=[];
	end %if
end %function

function fhat=get.fhat(h)
	if(h.fhat_is_set)
		fhat=nnfftmex('get_f_hat',h.plan);

		%switch h.d
		%case 1
			% Do nothing.
		%case 2
			% linearization in matlab with column (:) operator is columnwise, in NFFT it is rowwise
			%fhat=reshape(fhat,h.N_2,h.N_1).';
			%fhat=fhat(:);
		%case 3
			% linearization in matlab with column (:) operator is columnwise, in NFFT it is rowwise
			%fhat=reshape(fhat,h.N_3,h.N_2,h.N_1);
			%fhat=permute(fhat,[3,2,1]);
			%fhat=fhat(:);
		%otherwise
		%	error('Unknown error.');
		%end %switch
	else
		fhat=[];
	end %if
end %funcition

function f=get.f(h)
	if(h.f_is_set)
		f=nnfftmex('get_f',h.plan);
	else
		f=[];
	end %if
end %function

function N=get.N(h)
	N=[h.N_1;h.N_2;h.N_3];
end %function

function N=get.N1(h)
	N=[h.N1_1;h.N1_2;h.N1_3];
end %function

% User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nnfft_precompute_psi(h)
% Precomputations for NNFFT.
	if(~h.x_is_set)
		error('Before doing precomputations you have to set nodes in x.');
	else
		nnfftmex('precompute_psi',h.plan)
		h.precomputations_done=true;
	end %if
end %function

function nnfft_trafo_direct(h)
%% NNDFT.
%%
%% nnfft_trafo_direct(h)
%%
%% INPUT
%%   h  object of class type nnfft

	if(~h.precomputations_done)
		error('Before doing a NNFFT transform you have to do precomputations.');
	elseif(~h.fhat_is_set)
		error('Before doing a NNFFT transform you have to set Fourier coefficients in fhat.');
	else
		nnfftmex('trafo_direct',h.plan);
		h.f_is_set=true;
	end %if
end %function

function nnfft_trafo(h)
% NFFT.
%
% nfft_trafo(h)
%
% INPUT
%   h  object of class type nnfft

	if(~h.precomputations_done)
		error('Before doing a NNFFT transform you have to do precomputations.');
	elseif(~h.fhat_is_set)
		error('Before doing a NNFFT transform you have to set Fourier coefficients in fhat.');
	else
		nnfftmex('trafo',h.plan);
		h.f_is_set=true;
	end %if
end %function


end %methods

end %classdef

