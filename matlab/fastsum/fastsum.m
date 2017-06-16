% This class provides a Matlab interface to the NFFT-fastsum library.
%
% Examples
%   See Matlab scripts ...
classdef fastsum < handle

properties(Dependent=true)
	x {mustBeReal}        % source nodes (real Nxd matrix)
	alpha {mustBeNumeric} % fourier coefficients (column vector of length N)
	y {mustBeReal}        % target nodes (real Mxd matrix)
end %properties

properties(Dependent=true,SetAccess='private')
    f;      % function evaluations in y
    b;      % expansion coefficients
end %properties

properties(SetAccess='private')
	d {mustBeInteger,mustBePositive}    % spatial dimension 
    kernel
    c {mustBeNumeric}                   % kernel parameter
    flags {mustBeInteger,mustBeNonnegative} = 0
    n {mustBeInteger,mustBePositive}    % expansion degree (of NFFT)
    p {mustBeInteger,mustBePositive}    % degree of smoothness of regularization
    eps_I {mustBeNonnegative}           % inner boundary
    eps_B {mustBeNonnegative}           % outer boundary
    nn_x {mustBeInteger,mustBePositive} % oversampled nn in x
    nn_y {mustBeInteger,mustBePositive} % oversampled nn in y
    m_x {mustBeInteger,mustBePositive}  % NFFT-cutoff in x
    m_y {mustBeInteger,mustBePositive}  % NFFT-cutoff in y
end %properties

properties(Hidden=true,SetAccess='private',GetAccess='private')
	plan
	x_storage     % source nodes (real Nxd matrix)
	alpha_storage % fourier coefficients (column vector of length N)
    
	f_is_set=false             % flag if f is set
	plan_is_set=false          % flag if plan was created
	pre_x_done=false % flag if precomputations were done
	pre_y_done=false % flag if precomputations were done
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
    
    h.plan=fastsummex('init',h.d,h.kernel,h.c,h.flags,h.n,h.p,h.eps_I,h.eps_B);
    h.plan_is_set=true;
end %function

function delete(h)
% Destructor
	if(h.plan_is_set)
		fastsummex('finalize',h.plan);
	end %if
end %function

% Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set.x(h,x)
    h.x_storage = x;
    h.pre_x_done = false;
end %function

function set.alpha(h,alpha)
    h.alpha_storage = alpha;
    h.pre_x_done = false;
end %function

function set.y(h,y)
    fastsummex('set_y',h.plan,y,h.nn_y,h.m_y)
    h.pre_y_done = true;
end %function

function set.nn_x(h,nn_x)
    h.nn_x = nn_x + mod(nn_x,2);    % make nn even
end %function

function set.nn_y(h,nn_y)
    h.nn_y = nn_y + mod(nn_y,2);    % make nn even
end %function


% Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=get.x(h)
	x=h.x_storage;
end %function

function alpha=get.alpha(h)
	alpha=h.alpha_storage;
end %function

function y=get.y(h)
    if h.pre_y_done
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

% User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fastsum_trafo(h)
% NFFT.fastsum
%
% fastsum_trafo(h)
%
% INPUT
%   h  object of class type nfft
    if(~h.pre_x_done)
        fastsummex('set_x_alpha',h.plan,h.x,h.alpha,h.nn_x,h.m_x);
        h.pre_x_done=true;
    end
    
    if(~h.pre_y_done)
        fastsummex('set_y',h.plan,h.y,h.nn_y,h.m_y);
        h.pre_y_done=true;
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
    if(~h.pre_x_done)
        fastsummex('set_x_alpha',h.x,h.alpha,h.nn_x,h.m_x);
        h.pre_x_done=true;
    end
    
    if(~h.pre_y_done)
        fastsummex('set_y',h.plan,h.y,h.nn_y,h.m_y);
        h.pre_y_done=true;
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
        x.nn_x,s.m_x,s.nn_y,s.m_y);
    h.x=s.x;
    h.alpha=s.alpha;
    h.y=s.y;
end
end %methods
end %classdef