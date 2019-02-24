% This class provides a Matlab interface for the inverse adjoint NFFT.

classdef infft_adjoint < handle
    
    properties(Dependent=true)
        fhat                    % Given Fourier coefficients
    end 
    
    properties(SetAccess=private)
        y                       % Nodes, where inversion should be done
        N                       % Number of Fourier coefficients
        n                       % Expansion degree
        f                       % Approximated function values
        f_direct                % Directly computed function values
        m                       % Cut-off parameter
        sigma                   % Oversampling factor
        p                       % Degree of smoothness              (needed in quadratic setting M=N)
        eps_B                   % Outer boundary                    (needed in quadratic setting M=N)
        eps_I                   % Inner boundary (<=1/4)            (needed in quadratic setting M=N)
        c                       % Precomputed coefficients          (precomputed in quadratic setting M=N)
        d                       % Precomputed coefficients          (precomputed in quadratic setting M=N)
        m2                      % Cut-off parameter for inner NFFT  (needed in setting M~=N)
        window                  % Window function for inversion     (needed in setting M~=N)
        B_opt                   % Optimized matrix B*               (precomputed in setting M~=N)
        times                   % Computation times  
    end
    
    properties(GetAccess=private, SetAccess=private)
        M                       % Number of nodes
        y_storage               % Stored nodes (possibly sorted)
        fhat_storage            % Stored Fourier coefficients (possibly sorted)
        fhat_is_set = false;    % Flag if Fourier coefficients are set
        trafo_done = false;     % Flag if trafo is done
        direct_done = false;    % Flag if direct trafo is done
        perm                    % Permutation to sort y (if unsorted)
        x                       % Additional equidistant nodes      (needed in quadratic setting M=N)
        nn_oversampled          % oversampled n                     (needed in quadratic setting M=N)
        dist                    % Distance needed for a shift       (needed in quadratic setting M=N)
    end
    
    
    methods
        function h = infft_adjoint(y,N,varargin) % Constructor
           
           % Add further NFFT methods to search path
           s=fileparts(mfilename('fullpath'));
           s1=strcat(s,'/../fastsum');
           addpath(s1)
           s2=strcat(s,'/../nfft');
           addpath(s2)
            
           % Check input and determine whether nodes are in the correct interval
           if ( nargin==0 || sum(size(y))==0 )
               error('No nodes were given. The input "y" is compulsory.');
           elseif not(isvector(y))
               error('First argument has to be a vector.');
           elseif ( sum(y<-0.5) + sum(y>=0.5) ~= 0 )
               error('Nodes are in the wrong interval. Only nodes in [-0.5,0.5) are allowed.');
           elseif ( mod(nargin,2)==1 && isempty(varargin) )
               N = length(y); % Only nodes were given, so we assume quadratic setting. (for compatibility)
           elseif ( mod(nargin,2)==1 && isfloat(N) )
               error('Wrong number of arguments.')
           elseif ( mod(nargin,2)==1 && mod(length(varargin),2) )
               varargin = [N,varargin]; % Only nodes and optional input was given, so we assume quadratic setting. (for compatibility)
               N = length(y); 
           elseif ( sum(size(N))>2 || mod(N,2)~=0 || N<=0 || ischar(N) )
               error('Second argument has to be an even natural number.')
           end
            
           % Nodes have to be a row vector
           h.y = transpose(y(:));
           
           % Set number of nodes
           h.M = length(h.y);
           
           % Set number of Fourier coefficients
           h.N = N;
            
           % In quadratic setting nodes have to be sorted
           if (h.M==h.N && issorted(h.y)==0)
               [h.y_storage,h.perm] = sort(h.y);
           else
               h.y_storage = h.y;
               h.perm = [];
           end

           % Parse for optional input
           P = inputParser;
           P.KeepUnmatched = true;
           addParameter(P,'n',h.N, @(x) assert(mod(x,2)==0 && x>h.N,'Arguments have to be even natural numbers bigger than N.'))
           addParameter(P,'m',4, @(x) assert(mod(x,1)==0 && x>0,'Arguments must be natural numbers.'))
           addParameter(P,'sigma',2, @(x) assert(x>=1 && not(ischar(x)),'Arguments must be real numbers larger equal 1.'))
           
           if h.N == h.M
               addParameter(P,'p',4, @(x) assert(mod(x,1)==0 && x>=0 && x<=12,'Arguments have to be natural numbers smaller equal 12.'))
               addParameter(P,'eps_I',[], @(x) assert(x>0 && x<=1/4,'Only positive values smaller equal 1/4 are allowed.'))
               parse(P,varargin{:});

               % Set additional parameters
               h.m = P.Results.m;
               h.sigma = P.Results.sigma;   
               h.n = P.Results.n;
               h.p = P.Results.p;
               h.eps_I = P.Results.eps_I;
               if isempty(h.eps_I) % Set default value for eps_I if it is not set by the user
                   const = 4;
                   h.eps_I = const * h.p/h.n;
                   while (h.eps_I > 1/4)
                       const = const/2;
                       h.eps_I = const * h.p/h.n;
                   end
               end
               
               h.eps_B = 0;
               h.nn_oversampled = ceil(h.sigma*h.n);
               if mod(h.nn_oversampled,2)==1 % Make it even if necessary
                   h.nn_oversampled = h.nn_oversampled+1;
               end
           else
               addParameter(P,'m2',4, @(x) assert(mod(x,1)==0 && x>0,'Arguments must be natural numbers.'))
               addParameter(P,'window','Dirichlet', @(x) assert(iscategory(categorical({'BSpline' 'Gaussian' 'Sinc' 'Bessel' 'Dirichlet'}),x),'Possible values are ''BSpline'', ''Gaussian'', ''Sinc'', ''Bessel'' or ''Dirichlet''.'))
               parse(P,varargin{:});

               % Set additional parameters
               h.m = P.Results.m;
               h.sigma = P.Results.sigma;   
               h.n = P.Results.n;
               h.m2 = P.Results.m2;
               h.window = P.Results.window;
           end %if
           
           % Node-dependent precomputations
           if h.M == h.N
               precompute_quadratic(h)
           elseif h.M > h.N
               precompute_underdetermined(h)
           else
               precompute_overdetermined(h)
           end
           
        end %function
        
        
    % Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function set.fhat(h,fhat) % Set function values
            % Check whether the parameters match 
            if ( isvector(fhat) ~= 1 )
                error('Input fhat must be a vector.');
            elseif ( length(fhat) ~= h.N )
                error('The input vector fhat needs to be of length N.');
            end
            
            % If the nodes were sorted, also sort fhat
            if h.perm
                fhat = fhat(h.perm);
            end
            
            h.fhat_storage = fhat(:);
            h.fhat_is_set = true;
            h.trafo_done = false;
            h.direct_done = false;
        end %function
        
        
    % Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                     
        function fhat=get.fhat(h) % Get Fourier coefficients
            if not(h.fhat_is_set)
                error('No Fourier coefficients were set.')
            end
            
            % If y got sorted, use the inverse permutation to get back the order of the user
            if h.perm
                fhat = h.fhat_storage(h.perm);
            else 
                fhat = h.fhat_storage;
            end
        end %function
        
        function f=get.f(h) % Get approximations of the function values
            if(~h.trafo_done)
                error('No trafo was done.');
            else
                f = h.f;
            end
        end %function 
        
        function f_direct=get.f_direct(h) % Get directly computed function values
            if(~h.direct_done)
                error('No trafo was done.');
            else
                f_direct = h.f_direct;
            end
        end %function
        
        
    % User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        function infft_direct_adjoint(h)
        % Exact computation of an adjoint iNDFT, i.e., inversion of the adjoint nonequispaced Fourier matrix.
            if not(h.fhat_is_set)
                error('Trafo cannot be done. No Fourier coefficients were set.')
            end
            
            tic
            A = zeros(h.M,h.N);
            j = 1:h.M;
            for k = -h.N/2 : h.N/2-1
                A(:,k + h.N/2 + 1) = exp(2*pi*1i*k*h.y(j));
            end
            h.f_direct = A'\h.fhat;
            h.times.t_direct = toc;
            h.direct_done = true;
        end %function
        
        function infft_trafo_adjoint(h)
        % Fast computation of an adjoint iNFFT.   
            if not(h.fhat_is_set)
                error('Trafo cannot be done. No Fourier coefficients were set.')
            end
            
            if h.M == h.N
                trafo_quadratic(h)
            else
                trafo_rectangular(h)
            end
        end %function
    
    end %methods
    
    
    % Private methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Access=private)
        
        function precompute_quadratic(h)
        % Precomputations for the quadratic setting.
            % Suppress specific warnings (nodes are in an interval that is larger than normal)
            warning('off','fastsum:alphaDeleted')
            warning('off','fastsum:nodesOutsideBall')
            
            tic
            % Set constants for fastsum
            o = ones(h.M,1);
            % Set additional equispaced nodes
            h.x = -0.5:1/h.M:0.5-1/h.M;
            
            % Set vector for correction of sign and shift x if necessary
            [vec,h.dist] = infft.compute_sign(h.x,h.y_storage,h.M);
            h.x = h.x + h.dist/2;
            
          	% Computation of coefficients c_l 
            plan = fastsum(1,'log_sin',pi,0,h.n,h.p,h.eps_I,h.eps_B,h.nn_oversampled,h.m);
            c_abs = infft.compute_coeff(h.M,plan,h.x,h.y_storage,o);
            
            % Computation of coefficients d_j    
            d_abs = infft.compute_coeff(h.M,plan,h.y_storage,h.y_storage,o);
           	d_abs = -d_abs;
            
            % Stabilization
            s = min(abs(min(abs(d_abs))), abs(min(abs(c_abs))));
            c_abs = c_abs + s;
            d_abs = d_abs - s;
            c_abs = exp(c_abs);
            d_abs = exp(d_abs);
            
            % Perform correction of sign
            h.c = c_abs .* vec;
            h.d = d_abs .* repmat([-1;1],h.M/2,1);
            
            % Computation time
            h.times.t_precompute = toc;
        end %function
        
        function precompute_underdetermined(h)
        % Precomputations for the underdetermined setting.
            % Suppress specific warning (matrix may be ill-conditioned)
            warning('off','MATLAB:nearlySingularMatrix')
        
            % Initialize optimized matrix
            tic
            h.B_opt = sparse(h.M,h.n);

            for l = -h.n/2:h.n/2-1

              % Determine indices of non-zero entries
              a = 1:h.M;
              ind1 = ((-h.m+l)/h.n<=h.y_storage-1) & (h.y_storage-1<=(h.m+l)/h.n);
              ind2 = ((-h.m+l)/h.n<=h.y_storage) & (h.y_storage<=(h.m+l)/h.n);
              ind3 = ((-h.m+l)/h.n<=h.y_storage+1) & (h.y_storage+1<=(h.m+l)/h.n);
              a = [a(ind1),a(ind2),a(ind3)];

              % Computation of columns bl

              if not(isempty(a))
              L = transpose(-h.n/2:h.n/2-1);
              z = zeros(h.n,length(a));

              % Set nodes
              for b = 1:length(a)      
                  z(:,b) = L./h.n-h.y_storage(a(b));
              end

              switch h.window
                  case 'Dirichlet'
                      % Computation via Dirichlet kernel
                      L = 1/h.n .* sin((h.N-1)*pi*z)./sin(pi*z);
                  otherwise
                      % Computation via NFFT
                      h.y_storage = h.y_storage(:);
                      N_oversampled = ceil(h.sigma*h.N);
                      if mod(N_oversampled,2)==1 % If it's odd, make it even
                        N_oversampled = N_oversampled+1;
                      end

                      % Initialization
                      plan = nfft(1,h.N,h.n*length(a),N_oversampled,h.m2,bitor(PRE_PHI_HUT,PRE_PSI),FFTW_MEASURE); % Use of nfft_init_guru

                      % Set nodes
                      z = reshape(z,h.n*length(a),1);
                      z = mod(z+0.5,1)-0.5;
                      plan.x = -z; % NFFT defined with an additional minus

                      % Fourier coefficients
                      help = 1./(h.n*infft.phi_hat(h.window,h.m,h.N,h.n,-h.N/2:h.N/2-1));
                      plan.fhat = help(:);

                      % Node-dependent precomputation
                      nfft_precompute_psi(plan);

                      % Transformation
                      nfft_trafo(plan); 

                      % Get function values
                      L = plan.f;

                      % Reshape matrix
                      L = reshape(L,h.n,length(a));
              end

              % Solve normal equations
              vec = (L'*L) \ transpose(L(l+h.n/2+1,:));

              % Compose the whole matrix
              h.B_opt = h.B_opt + sparse(a,l+h.n/2+1,vec,h.M,h.n);
              end
            end
            
            h.B_opt = h.B_opt';
            h.times.t_precompute = toc;
        end %function
        
        function precompute_overdetermined(h)
        % Precomputations for the overdetermined setting.
            % Suppress specific warning (matrix may be ill-conditioned)
            warning('off','MATLAB:nearlySingularMatrix')
        
            % Initialize optimized matrix
            tic
            h.B_opt = sparse(h.n,h.M);

            for j = 1:h.M

              % Determine indices of non-zero entries
              a = mod((ceil(h.n*h.y_storage(j))-h.m : floor(h.n*h.y_storage(j))+h.m) + h.n/2 + 1, h.n);
              ind = (a==0);    
              a(ind) = h.n;       % 0 equivalent m modulo m

              % Computation of columns bj
              if not(isempty(a))
              L = -h.n/2:h.n/2-1;
              z = zeros(h.M,length(a));

              % Set nodes
              for b = 1:length(a)      
                  z(:,b) = h.y_storage-L(a(b))/h.n;
              end

              switch h.window
                  case 'Dirichlet'
                      % Computation via Dirichlet kernel
                      K = 1/h.n .* sin((h.N-1)*pi*z)./sin(pi*z);
                  otherwise
                      % Computation via NFFT
                      N_oversampled = ceil(h.sigma*h.N);
                      if mod(N_oversampled,2)==1 % If it's odd, make it even
                        N_oversampled = N_oversampled+1;
                      end

                      % Initialization
                      plan = nfft(1,h.N,length(a)*h.M,N_oversampled,h.m2,bitor(PRE_PHI_HUT,PRE_PSI),FFTW_MEASURE); % Use of nfft_init_guru

                      % Set nodes
                      z = reshape(z,h.M*length(a),1);
                      z = mod(z+0.5,1)-0.5;
                      plan.x = -z; % NFFT defined with an additional minus

                      % Fourier coefficients
                      help = 1./(h.n*infft.phi_hat(h.window,h.m,h.N,h.n,-(-h.N/2:h.N/2-1)));
                      plan.fhat = help(:);

                      % Node-dependent precomputation
                      nfft_precompute_psi(plan);

                      % Transformation
                      nfft_trafo(plan); 

                      % Get function values
                      K = plan.f;

                      % Reshape matrix
                      K = reshape(K,h.M,length(a));
              end

              % Solve normal equations
              vec = (K'*K) \ (h.N*transpose(K(j,:)));

              % Compose the whole matrix
              h.B_opt = h.B_opt + sparse(a,j,transpose(vec),h.n,h.M);
              end  
            end
            
            % Computation time
            h.times.t_precompute = toc;
        end %function
        
        function trafo_quadratic(h)
        % Computation of an adjoint iNFFT for the quadratic setting.             
            % Suppress specific warnings (nodes are in an interval that is larger than normal)
            warning('off','fastsum:alphaDeleted')
            warning('off','fastsum:nodesOutsideBall')
             
            tic
            % Computation of the Fourier coefficients by iFFT
            v = ifft(h.fhat_storage);
            v = v.*transpose(exp(pi*1i*(0:h.M-1)));
            v = fftshift(v);
            v = v.*transpose(exp(pi*1i*(-h.M/2:h.M/2-1)*h.dist));
                        
            % Set constants for fastsum of f
            alpha = v.*h.c;
            
            % Computation of coefficients f
            plan = fastsum(1,'cot',-pi,0,h.n,h.p,h.eps_I,h.eps_B,h.nn_oversampled,h.m);
            ftilde = infft.compute_coeff(h.M,plan,h.y_storage,h.x,alpha);
            
            % Final computation
            ftilde = h.d.*(ftilde + 1i*sum(alpha));
            
            % Set flag
            h.trafo_done = true;
            
            % If y got sorted, use the inverse permutation to get back the order of the user
            if h.perm
                h.f(h.perm) = ftilde;
            else 
                h.f = ftilde;
            end
            
            % Computation time
            h.times.t_trafo = toc;
        end %function
        
        function trafo_rectangular(h)
        % Computation of an adjoint iNFFT for the rectangular setting.
            tic
            % Perform an NFFT
            % Multiplication with diagonal matrix
            ftilde = transpose(1./infft.phi_hat(h.window,h.m,h.N,h.n,-h.N/2:h.N/2-1)) .* h.fhat;
            
            % Perform an iFFT
            ftilde = ifft([zeros((h.n-h.N)/2,1);ftilde;zeros((h.n-h.N)/2,1)]);
            ftilde = ftilde.*transpose(exp(pi*1i*(0:h.n-1)));
            ftilde = fftshift(ftilde);
            
            % Multiplication with adjoint of optimized sparse matrix
            ftilde = h.B_opt'*ftilde;
            
            % Multiplication with specific constant
            if h.M > h.N
                const = 1;
            else 
                const = 1/h.N;
            end
            ftilde = const.*ftilde;

            % Set approximations
            h.f = ftilde;

            % Computation time
            h.times.t_trafo = toc;
            
            % Set flag
            h.trafo_done = true;
        end %function
    end %methods
end %classdef