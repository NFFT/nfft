classdef infft < handle
    
    properties(Dependent=true)
        f                       % Dummy variable
    end
    
    properties(SetAccess=private)
        fbar                    % Approximated Fourier coefficients
        fbar_direct             % Direcly computed Fourier coefficients
        times                   % Computing times
    end
    
    properties(GetAccess=private, SetAccess=private)
        N                       % Number of knots
        x                       % Additional equidistant knots
        y                       % Knots where inversion should be done
        perm                    % Permutaion to sort y (if unsorted)
        f_storage               % Evaluations at points y
        abs                     % Distance needed for a shift
        trafo_done = false;     % Flag if trafo is done
        direct_done = false;    % Flag if direct trafo is done
        m                       % Cut-off parameter for NFFT 
        p                       % Degree of smoothness of fastsum
        sigma                   % Oversampling factor for NFFT
        n                       % Expansion degree  
        nn_oversampled          % n*sigma
        eps_B = 0;              % Outer boundary
        eps_I                   % Inner boundary (<=1/4)
        c                       % Needed coefficients for iNFFT
        d                       % Needed coefficients for iNFFT
    end
    
    methods
        function h = infft(y,varargin)
           % Constructor
           
           s=fileparts(mfilename('fullpath'));
           s=strcat(s,'/../fastsum');
           addpath(s)         
            
           % Check input and determine whether knots are in the correct interval
           if ( nargin==0 || sum(size(y))==0 )
               error('No knots were given. The input "y" is compulsory.');
           elseif mod(nargin,2)==0
               error('Wrong number of arguments.')
           elseif ( size(y,1)~=1 && size(y,2)~=1 )
               error('Input must be a vector.');
           elseif ( sum(y < -0.5) + sum(y >= 0.5) ~= 0 )
               error('Knots are in the wrong interval. Only knots in [-0.5,0.5) are allowed.');
           end
            
           % Make it a row vector if necessary
           if size(y,2)==1             
               h.y = transpose(y);
           else
               h.y = y;
           end 
            
           % Sort it if necessary
           if issorted(h.y)==0
               [h.y,h.perm] = sort(h.y);
           else
               h.perm = [];
           end
           
           % Number of knots
           h.N = length(h.y);

           % Optional Input
           P = inputParser;
           P.KeepUnmatched = true;
           
           addParameter(P,'m',4, @(x) (x>0 && round(x)==x && isreal(x)))
           addParameter(P,'p',4, @(x) (x>=0 && x<=12 && isreal(x)))
           addParameter(P,'sigma',2, @(x) (x>=1 && isreal(x)))
           addParameter(P,'n',2*h.N, @(x) (x>0 && mod(x,2) == 0 && isreal(x)))
           addParameter(P,'eps_I',NaN, @(x) (x>0 && x<=1/4 && isreal(x)))
                    
           parse(P,varargin{:});

           h.m = P.Results.m;
           h.p = P.Results.p;
           h.sigma = P.Results.sigma;   
           h.n = P.Results.n;
           h.eps_I = P.Results.eps_I;
           
           h.nn_oversampled = ceil(h.sigma*h.n);             
           % Check if it is odd and make it even if necessary            
           if mod(h.nn_oversampled,2)==1
               h.nn_oversampled = h.nn_oversampled+1;
           end
           
           % Set default value for eps_I if it is not set
           if not(isfinite(h.eps_I))
               const = 2;
               h.eps_I = const * 2*h.p/h.n;
               while (h.eps_I > 1/4)
                   const = const/2;
                   h.eps_I = const * 2*h.p/h.n;
               end
           end
           
           % Suppress warnings (knots are in an interval that is bigger than normal)
           warning('off','fastsum:alphaDeleted')
           warning('off','fastsum:nodesOutsideBall')
           
           tic
           % Set constants for fastsum
           o = ones(h.N,1);
           % Set additional equidistant knots
           h.x = -0.5:1/h.N:0.5-1/h.N;
           
           % Set vector for correction of sign and shift x if necessary
           vec = compute_sign(h);
           help = toc;
           
           % Computation of coefficients c_l 
           tic
           plan = fastsum(1,'log_sin',pi,0,h.n,h.p,h.eps_I,h.eps_B,h.nn_oversampled,h.m);
           h.c = compute_coeff(h,plan,h.x,h.y,o);
           h.times.t_c = toc;

           % Computation of coefficients d_j    
           tic
           h.d = compute_coeff(h,plan,h.y,h.y,o);
           h.d = -h.d;
           h.times.t_d = toc;
           
           % Set constant and divide it out
           tic
           s = min(abs(min(abs(h.d))), abs(min(abs(h.c))));
           h.c = h.c + s;
           h.d = h.d - s;
           h.c = exp(h.c);
           h.d = exp(h.d);
           
           % Correction of sign
           h.c = h.c .* vec;
           h.d = h.d .* repmat([-1;1],h.N/2,1);
           
           h.times.t_total = help + toc;
        end %function
        
    % Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function set.f(h,f)
            % Check whether the parameters y and f match 
            if ( isvector(f) ~= 1 )
                error('Input must be a vector.');
            elseif ( length(f) ~= h.N )
                error('The input vectors y and f need to have same length.');
            end
            
            % Make it a column vector
            f = f(:);
            
            % If y needed to get sorted, also sort f
            if h.perm
                f = f(h.perm);
            end
            
            h.f_storage = f;
            h.trafo_done = false;
            h.direct_done = false;
        end %function
        
    % Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                     
        function fbar=get.fbar(h)
            if(~h.trafo_done)
                error('No trafo was done.');
            else
                fbar = h.fbar;
            end
        end %function 
        
        function fbar_direct=get.fbar_direct(h)
            if(~h.direct_done)
                error('No trafo was done.');
            else
                fbar_direct = h.fbar_direct;
            end
        end %function
        
    % User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        function infft_direct(h,f)
        % Exact computation of an iNFFT, i.e. inversion of the nonequispaced Fourier-matrix.
            tic
            A = zeros(h.N,h.N);
            j = 1:h.N;
            for k = -h.N/2 : h.N/2-1
                A(:,k + h.N/2 + 1) = exp(2*pi*1i*k*h.y(j));
            end
            h.fbar_direct = A\f;
            h.times.t_direct = toc;
            h.direct_done = true;
        end %function
        
        function infft_trafo(h)
        % Final computations for iNFFT.             
            % Set constants for fastsum of g_l
            tic
            alpha = h.f_storage.*h.d;
            h.times.t_total = h.times.t_total + toc;
            
            % Computation of coefficients g_l
            tic
            plan = fastsum(1,'cot',pi,0,h.n,h.p,h.eps_I,h.eps_B,h.nn_oversampled,h.m);
            g = compute_coeff(h,plan,h.x,h.y,alpha);
            h.times.t_g = toc; 
            
            % Final computation
            tic
            g = -h.c.*(-g + 1i*h.d'*h.f_storage);
            h.times.t_total = h.times.t_total + toc;
            
            % Computation of the Fourier-coefficients by FFT
            tic
            fhat = fft(g);
            fhat = fhat.*transpose(exp(-pi*1i*(0:h.N-1)));
            fhat = fftshift(fhat);
            fhat = fhat/h.N;   
            fhat = fhat.*transpose(exp(-pi*1i*(-h.N/2:h.N/2-1)*h.abs));
            h.times.t_fft = toc;
            
            % Time required
            h.times.t_total = h.times.t_total + h.times.t_c + h.times.t_d + h.times.t_g + h.times.t_fft;
                        
            % If y got sorted, use the inverse permutation to get back the order of the user
            if h.perm
                h.fbar(h.perm) = fhat;
            else 
                h.fbar = fhat;
            end
            
            h.trafo_done = true;
        end %function
    
    end %methods
    
    % Private methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Access=private)
        function vec = compute_sign(h)
            % Auxiliary function for computing the correction of sign.
            anz = zeros(h.N,1);
            h.abs = 1;
            i = 1;
            j = 1;

            while (i <= h.N && j <= h.N)
            if h.x(i) < h.y(j)
                diff = h.y(j)-h.x(i);
                if diff < h.abs
                    h.abs = diff;
                end
                i = i + 1;
                if i <= h.N
                    anz(i) = anz(i-1);
                end
            elseif h.x(i) > h.y(j)
                anz(i) = anz(i) + 1;
                j = j + 1;
            else
                if j < h.N
                    diff = h.y(j+1)-h.x(i);
                    if diff < h.abs
                        h.abs = diff;
                    end
                end
                i = i + 1;
                if i <= h.N
                    anz(i) = anz(i-1);
                end
                anz(i-1) = anz(i-1) + 1;
            end
            end

            % If necessary do a shift
            if not(isempty(intersect(h.x,h.y)))
                h.x = h.x + h.abs/2;
            else
                h.abs = 0;
            end

            vec = (-1).^anz;
        end %function
        
        function coeff = compute_coeff(h,plan,a,b,alpha)
            % Auxiliary function for computing coefficients.
            coeff = zeros(h.N,1);

            % New index sets for knots a
            A1 = a<-0.25;
            A2 = a>=0.25;
            A3 = and(not(A1),not(A2));

            % New index sets for knots b
            B1 = b<0.25;
            B2 = not(B1);
            B3 = b>=-0.25;
            B4 = not(B3);

            %------------------------------------------------------------------
            % Case 1 (a < -0.25)
            % Without shift of knots

            plan.x = b(B1)';
            plan.alpha = alpha(B1);
            plan.y = a(A1)';
            fastsum_trafo(plan)
            coeff(A1) = plan.f;

            % With shift

            plan.x = (b(B2)-1/4)';
            plan.alpha = alpha(B2);
            plan.y = (a(A1)+3/4)';
            fastsum_trafo(plan)
            coeff(A1) = coeff(A1) + plan.f;
            
            %------------------------------------------------------------------
            % Case 2 (a >= 0.25)
            % Without shift of knots

            plan.x = b(B3)';
            plan.alpha = alpha(B3);
            plan.y = a(A2)';
            fastsum_trafo(plan)
            coeff(A2) = plan.f;

            % With shift

            plan.x = (b(B4)+1/4)';
            plan.alpha = alpha(B4);
            plan.y = (a(A2)-3/4)';
            fastsum_trafo(plan)
            coeff(A2) = coeff(A2) + plan.f;

            %------------------------------------------------------------------
            % Case 3 (-0.25 <= a < 0.25)

            plan.x = b';
            plan.alpha = alpha;
            plan.y = a(A3)';
            fastsum_trafo(plan)
            coeff(A3) = plan.f;
            
        end %function
    end % methods
end %classdef