
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

% This class provides a Matlab interface to the FASTSUM application.
classdef fastsum < handle

properties(Dependent=true)
	x           % source nodes in d-ball with radius 1/4-eps_b/2 (real Nxd matrix)
	alpha       % fourier coefficients (column vector of length N)
	y           % target nodes in d-ball with radius 1/4-eps_b/2 (real Mxd matrix)
end %properties

properties(Dependent=true,SetAccess='private')
    f;          % function evaluations in y
    b;          % expansion coefficients
    N;          % number of source nodes
end %properties

properties(SetAccess='private')
	d               % spatial dimension 
    kernel          % kernel
    c               % kernel parameter
    flags  = 0      % flags
    n               % expansion degree (of NFFT)
    p               % degree of smoothness of regularization
    eps_I           % inner boundary
    eps_B           % outer boundary
    nn_x            % oversampled nn in x
    nn_y            % oversampled nn in y
    m_x             % NFFT-cutoff in x
    m_y             % NFFT-cutoff in y
end %properties

properties(Hidden=true,SetAccess='private',GetAccess='private')
	plan
	x_is_set=false      % flag if x is set
	y_is_set=false      % flag if y is set
	alpha_is_set=false  % flag if alpha is set
	f_is_set=false      % flag if f is set
	plan_is_set=false   % flag if plan was created
end %properties

methods

% Constructer and destructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=fastsum(d,kernel,c,flags,n,p,eps_I,eps_B,nn_x,m_x,nn_y,m_y)
% Constructor
%
% h=fastsum(d,kernel,c,flags,n,p,eps_I,eps_B)
% h=fastsum(d,kernel,c,flags,n,p,eps_I,eps_B,nn,m)
% h=fastsum(d,kernel,c,flags,n,p,eps_I,eps_B,nn_x,m_x,nn_y,m_y)
%
% INPUT
% kernel = 'multiquadric', etc. (see options below)
% c       kernel parameter
% n       expansion degree
% m       cut-off parameter for NFFT
% p       degree of smoothness of regularization
% eps_I   inner boundary
% eps_B   outer boundary
% nn      oversampled nn in NFFT (in x and y)
% m       cut-off parameter for NFFT (in x and y)
% nn_x    oversampled nn in NFFT (in x)
% m_x     cut-off parameter for NFFT (in x)
% nn_y    oversampled nn in NFFT (in y)
% m_y     cut-off parameter for NFFT (in y)
%
% OUTPUT
%   h   object of class type fastsum
    narginchk(8,12);

	h.d = d;
    h.kernel = kernel;
    h.c = c;
    h.flags = flags;
    h.n = n;
    h.p = p;
    h.eps_I = eps_I;
    h.eps_B = eps_B;
    
    switch nargin
    case 8
        h.nn_x = 2*n;   % default oversampling factor 2
        h.nn_y = 2*n;
        h.m_x = p;      % default NFFT-cutoff
        h.m_y = p;
    case 10
        h.nn_x=nn_x;
        h.nn_y=nn_x;
        h.m_x=m_x;
        h.m_y=m_x;
    case 12
        h.nn_x=nn_x;
        h.nn_y=nn_y;
        h.m_x=m_x;
        h.m_y=m_y;
    otherwise
        error('Wrong number of arguments')
    end
    
    if nn_x<=h.n || nn_x<=h.n
        error('Oversampled nn must be larger than n')
    end
    
    h.plan=fastsummex('init',h.d,h.kernel,h.c,h.flags,h.n,h.p,h.eps_I,h.eps_B);
    h.plan_is_set=true;
end %function

function delete(h)
% Destructor
	if(h.plan_is_set)
		fastsummex('finalize',h.plan);
		h.plan_is_set=false;
	end %if
end %function

% Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set.x(h,x)
    if (size(x,2)~=h.d || ~isreal(x))
        error('x must be a real matrix with d columns');
    end
    if h.x_is_set && h.alpha_is_set && size(x,1)~=h.N
        warning('fastsum:alphaDeleted',...
            'Size of x has changed, so alpha is deleted');
        h.alpha_is_set=false;
    end
    fastsummex('set_x',h.plan,x,h.nn_x,h.m_x)
    h.x_is_set = true;
end %function

function set.alpha(h,alpha)
    if ~h.x_is_set
        error('You have to set x before you set alpha')
    end
    if (~iscolumn(alpha) || ~isnumeric(alpha) || length(alpha)~=h.N)
        error('alpha must be a column vector of length N');
    end
    fastsummex('set_alpha',h.plan,alpha)
    h.alpha_is_set=true;
end %function

function set.y(h,y)
    if (size(y,2)~=h.d || ~isreal(y))
        error('y must be a real matrix with d columns');
    end
    fastsummex('set_y',h.plan,y,h.nn_y,h.m_y)
    h.y_is_set = true;
end %function

function set.d(h,d)
    if (~isscalar(d) || ~isreal(d) || d<=0 || mod(d,1)~=0)
        error('d must be a positive integer')
    end
    h.d = d;
end %function

function set.flags(h,n)
    if (~isscalar(n) || ~isreal(n) || n<0 || mod(n,1)~=0)
        error('flags must be a non-negative integer')
    end
    h.flags = n;
end %function

function set.n(h,n)
    if (~isscalar(n) || ~isreal(n) || n<=0 || mod(n,2)~=0)
        error('n must be a positive, even integer')
    end
    h.n = n;
end %function

function set.p(h,p)
    if (~isscalar(p) || ~isreal(p) || p<=0 || mod(p,1)~=0)
        error('p must be a positive integer')
    end
    h.p = p;
end %function

function set.eps_I(h,arg)
    if (~isscalar(arg) || arg<0 || arg>=pi/4)
        error('eps_I must be in the interval [0,pi/4)')
    end
    h.eps_I = arg;
end %function

function set.eps_B(h,arg)
    if (~isscalar(arg) || arg<0 || arg>=pi/4)
        error('eps_B must be in the interval [0,pi/4)')
    end
    h.eps_B = arg;
end %function

function set.nn_x(h,arg)
    if (~isscalar(arg) || ~isreal(arg) || arg<=0 || mod(arg,1)~=0)
        error('nn_x must be a positive integer')
    end
    h.nn_x = arg;
end %function

function set.m_x(h,arg)
    if (~isscalar(arg) || ~isreal(arg) || arg<=0 || mod(arg,1)~=0)
        error('m_x must be a positive integer')
    end
    h.m_x = arg;
end %function

function set.nn_y(h,arg)
    if (~isscalar(arg) || ~isreal(arg) || arg<=0 || mod(arg,1)~=0)
        error('nn_y must be a positive integer')
    end
    h.nn_y = arg;
end %function

function set.m_y(h,arg)
    if (~isscalar(arg) || ~isreal(arg) || arg<=0 || mod(arg,1)~=0)
        error('m_y must be a positive integer')
    end
    h.m_y = arg;
end %function


% Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=get.x(h)
    if h.x_is_set
        x=fastsummex('get_x',h.plan);
    else
        x=[];
    end
end %function

function alpha=get.alpha(h)
    if h.alpha_is_set
        alpha=fastsummex('get_alpha',h.plan);
    else
        alpha=[];
    end
end %function

function y=get.y(h)
    if h.y_is_set
        y=fastsummex('get_y',h.plan);
    else
        y=[];
    end
end %function

function b=get.b(h)
	if(h.plan_is_set)
		b=fastsummex('get_b',h.plan);
	else
		b=[];
	end %if
end %function

function f=get.f(h)
	if(h.f_is_set)
		f=fastsummex('get_f',h.plan);
	else
		f=[];
	end %if
end %function

function N=get.N(h)
	if(h.x_is_set)
		N=fastsummex('get_N_total',h.plan);
	else
		N=[];
	end %if
end %function

% User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fastsum_trafo(h)
% NFFT.fastsum
%
% fastsum_trafo(h)
%
% INPUT
%   h  object of class type nfft
    if(~h.x_is_set)
		error('Before doing a fastsum transform you have to set nodes x.');
    end
    if(~h.alpha_is_set)
		error('Before doing a fastsum transform you have to set coefficients alpha.');
    end
    if(~h.y_is_set)
		error('Before doing a fastsum transform you have to set nodes y.');
    end
        
	fastsummex('trafo',h.plan);
	h.f_is_set=true;
end %function

function fastsum_trafo_direct(h)
% NFFT.fastsum
%
% fastsum_trafo_direct(h)
%
% INPUT
%   h  object of class type nfft
    if(~h.x_is_set)
		error('Before doing a fastsum transform you have to set nodes x.');
    end
    if(~h.alpha_is_set)
		error('Before doing a fastsum transform you have to set coefficients alpha.');
    end
    if(~h.y_is_set)
		error('Before doing a fastsum transform you have to set nodes y.');
    end
    
    fastsummex('trafo_direct',h.plan);
    h.f_is_set=true;
end %function

% Save and load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = saveobj(h)
    s=struct('x',h.x,'y',h.y,'alpha',h.alpha,'d',h.d,'kernel',h.kernel,'c',h.c,...
        'flags',h.flags,'n',h.n,'p',h.p,'eps_I',h.eps_I,'eps_B',h.eps_B,...
        'nn_x',h.nn_x,'nn_y',h.nn_y,'m_x',h.m_x,'m_y',h.m_y);
end

end %methods

methods(Static)
function h = loadobj(s)
    h=fastsum(s.d,s.kernel,s.c,s.flags,s.n,s.p,s.eps_I,s.eps_B,...
        s.nn_x,s.m_x,s.nn_y,s.m_y);
    if ~isempty(s.x)
        h.x=s.x;
    end
    if ~isempty(s.alpha)
        h.alpha=s.alpha;
    end
    if ~isempty(s.y)
        h.y=s.y;
    end
end
end %methods
end %classdef
