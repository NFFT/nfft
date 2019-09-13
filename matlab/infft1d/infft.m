% This class provides a Matlab interface for the iNFFT.

classdef infft < handle
    
    properties(Dependent=true)
        f                       % Given function values
        fhat                    % Given Fourier coefficients
    end 
    
    properties(SetAccess=private)
        y                       % Nodes, where inversion should be done
        N                       % Number of Fourier coefficients
        n                       % Expansion degree
        fcheck                  % Approximated Fourier coefficients
        fcheck_direct           % Directly computed Fourier coefficients
        ftilde                  % Approximated function values
        ftilde_direct           % Directly computed function values
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
        flag_toeplitz           % Flag for exact computation via Toeplitz approach (only possible for M>N)
        D                       % Components of Gohberg-Semencul formula (precomputed if flag_toeplitz==1)
        times                   % Computation times  
    end
    
    properties(GetAccess=private, SetAccess=private)
        M                               % Number of nodes
        y_storage                       % Stored nodes (possibly sorted)
        f_storage                       % Stored function values (possibly sorted)
        fhat_storage                    % Stored Fourier coefficients (possibly sorted)
        f_is_set = false;               % Flag if function values are set
        fhat_is_set = false;            % Flag if Fourier coefficients are set
        trafo_done = false;             % Flag if trafo is done
        direct_done = false;            % Flag if direct trafo is done
        trafo_adjoint_done = false;     % Flag if trafo is done
        direct_adjoint_done = false;    % Flag if direct trafo is done
        perm                            % Permutation to sort y (if unsorted)
        x                               % Additional equidistant nodes      (needed in quadratic setting M=N)
        nn_oversampled                  % oversampled n                     (needed in quadratic setting M=N)
        dist                            % Distance needed for a shift       (needed in quadratic setting M=N)
    end
    
    
    methods
        function h = infft(y,N,varargin) % Constructor
           
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
           addParameter(P,'flag_toeplitz',0, @(x) assert((x==0 || x==1) && h.M>=h.N,'Arguments must be zero or one. Besides, computation via Gohberg-Semencul formula is only possible in case M>=N.'))
           
           if h.N == h.M
               addParameter(P,'p',4, @(x) assert(mod(x,1)==0 && x>=0 && x<=12,'Arguments have to be natural numbers smaller equal 12.'))
               addParameter(P,'eps_I',[], @(x) assert(x>0 && x<=1/4,'Only positive values smaller equal 1/4 are allowed.'))
               parse(P,varargin{:});

               % Set additional parameters
               h.m = P.Results.m;
               h.sigma = P.Results.sigma;
               h.flag_toeplitz = P.Results.flag_toeplitz;
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
               addParameter(P,'window','Dirichlet', @(x) assert(iscategory(categorical({'BSpline' 'Gaussian' 'Sinc' 'Kaiser' 'Bessel' 'Dirichlet'}),x),'Possible values are ''BSpline'', ''Gaussian'', ''Sinc'', ''Kaiser'', ''Bessel'' or ''Dirichlet''.'))
               parse(P,varargin{:});

               % Set additional parameters
               h.m = P.Results.m;
               h.sigma = P.Results.sigma;   
               h.n = P.Results.n;
               h.m2 = P.Results.m2;
               h.window = P.Results.window;
               h.flag_toeplitz = P.Results.flag_toeplitz;
           end %if
           
           % Node-dependent precomputations
           if h.flag_toeplitz == 1
               precompute_toep(h)
           elseif h.M == h.N
               precompute_quadratic(h)
           elseif h.M > h.N
               precompute_overdetermined(h)
           else
               precompute_underdetermined(h)
           end
           
        end %function
        
        
    % Set functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function set.f(h,f) % Set function values
            % Check whether the parameters y and f match 
            if ( isvector(f) ~= 1 )
                error('Input f must be a vector.');
            elseif ( length(f) ~= h.M )
                error('The input vectors y and f need to be of same length.');
            end
            
            % If the nodes were sorted, also sort f
            if h.perm
                f = f(h.perm);
            end
            
            h.f_storage = f(:);
            h.f_is_set = true;
            h.trafo_done = false;
            h.direct_done = false;
        end %function
        
        function set.fhat(h,fhat) % Set Fourier coefficients
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
            h.trafo_adjoint_done = false;
            h.direct_adjoint_done = false;
        end %function
        
        
    % Get functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                     
        function f=get.f(h) % Get function values
            if not(h.f_is_set)
                error('No function values were set.')
            end
            
            % If y got sorted, use the inverse permutation to get back the order of the user
            if h.perm
                f = h.f_storage(h.perm);
            else 
                f = h.f_storage;
            end
        end %function
        
        function fcheck=get.fcheck(h) % Get approximations of the Fourier coefficients
            if(~h.trafo_done)
                error('No trafo was done.');
            else
                fcheck = h.fcheck;
            end
        end %function 
        
        function fcheck_direct=get.fcheck_direct(h) % Get directly computed coefficients
            if(~h.direct_done)
                error('No trafo was done.');
            else
                fcheck_direct = h.fcheck_direct;
            end
        end %function
        
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
        
        function ftilde=get.ftilde(h) % Get approximations of the function values
            if(~h.trafo_done)
                error('No trafo was done.');
            else
                ftilde = h.ftilde;
            end
        end %function 
        
        function ftilde_direct=get.ftilde_direct(h) % Get directly computed function values
            if(~h.direct_done)
                error('No trafo was done.');
            else
                ftilde_direct = h.ftilde_direct;
            end
        end %function

        
    % User methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        function infft_trafo(h)
        % Fast computation of an iNFFT.   
            if not(h.f_is_set)
                error('Trafo cannot be done. No function values were set.')
            end
            
            if h.flag_toeplitz == 1
                trafo_toep(h)
            elseif h.M == h.N
                trafo_quadratic(h)
            else
                trafo_rectangular(h)
            end
        end %function
        
        function infft_trafo_direct(h)
        % Exact computation of an iNDFT, i.e., inversion of the nonequispaced Fourier matrix.
            if not(h.f_is_set)
                error('Trafo cannot be done. No function values were set.')
            end
            
            tic
            A = zeros(h.M,h.N);
            j = 1:h.M;
            for k = -h.N/2 : h.N/2-1
                A(:,k + h.N/2 + 1) = exp(2*pi*1i*k*h.y(j));
            end
            h.fcheck_direct = A\h.f;
            h.times.t_direct = toc;
            h.direct_done = true;
        end %function
        
        function infft_direct(h)
        % Additional function for compatibility
            infft_trafo_direct(h)
        end %function
        
        function infft_adjoint(h)
        % Fast computation of an adjoint iNFFT.   
            if not(h.fhat_is_set)
                error('Trafo cannot be done. No Fourier coefficients were set.')
            end
            
            if h.flag_toeplitz == 1
                error('Trafo cannot be done. Computation via Gohberg-Semencul formula is not available in the adjoint setting.')
            end
            
            if h.M == h.N
                adjoint_quadratic(h)
            else
                adjoint_rectangular(h)
            end
        end %function
        
        function infft_adjoint_direct(h)
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
            h.ftilde_direct = A'\h.fhat;
            h.times.t_direct = toc;
            h.direct_done = true;
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
        
        function precompute_overdetermined(h)
        % Precomputations for the overdetermined setting.
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
        
        function precompute_underdetermined(h)
        % Precomputations for the underdetermined setting.
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
        
        function precompute_toep(h)
        % Precomputations in case flag_toeplitz is true.
            tic
            % Initialize matrix A
            A = zeros(h.M,h.N);
            for k = -h.N/2 : h.N/2-1
                A(:,k + h.N/2 + 1) = exp(2*pi*1i*k*h.y);
            end

            % Compute inverse of Toeplitz matrix A'*A
            [M1_c,M1_r,M2_c,M2_r,M3_c,M3_r,M4_c,M4_r] = infft.ToeplitzInv(A'*A(:,1),(A'*A(:,1))');
            
            % Compute diagonalization of matrices M1-M4
            h.D.D1 = infft.ToeplitzDiag(M1_c,M1_r);
            h.D.D2 = infft.ToeplitzDiag(M2_c,M2_r);
            h.D.D3 = infft.ToeplitzDiag(M3_c,M3_r);
            h.D.D4 = infft.ToeplitzDiag(M4_c,M4_r);
            h.D.x0 = M1_c(1);
            
            % Computation time
            h.times.t_precompute = toc;
        end %function
        
        function trafo_quadratic(h)
        % Computation of an iNFFT for the quadratic setting.             
            % Suppress specific warnings (nodes are in an interval that is larger than normal)
            warning('off','fastsum:alphaDeleted')
            warning('off','fastsum:nodesOutsideBall')
             
            tic
            % Set constants for fastsum of g
            alpha = h.f_storage.*h.d;
            
            % Computation of coefficients g
            plan = fastsum(1,'cot',pi,0,h.n,h.p,h.eps_I,h.eps_B,h.nn_oversampled,h.m);
            g = infft.compute_coeff(h.M,plan,h.x,h.y_storage,alpha);
            
            % Final computation
            g = -h.c.*(-g + 1i*sum(alpha));
            
            % Computation of the Fourier coefficients by FFT
            temp = fft(g);
            temp = temp.*transpose(exp(-pi*1i*(0:h.M-1)));
            temp = fftshift(temp);
            temp = temp/h.M;   
            temp = temp.*transpose(exp(-pi*1i*(-h.M/2:h.M/2-1)*h.dist));
            
            % Set flag
            h.trafo_done = true;
                        
            % If y got sorted, use the inverse permutation to get back the order of the user
            if h.perm
                h.fcheck(h.perm) = temp;
            else 
                h.fcheck = temp;
            end
            
            % Computation time
            h.times.t_trafo = toc;
        end %function
        
        function trafo_rectangular(h)
        % Computation of an iNFFT for the rectangular setting.
            tic
            % Perform an adjoint NFFT
            % Multiplication with optimized sparse matrix
            temp = h.B_opt*h.f_storage;

            % Perform an FFT
            temp = fft(temp);
            temp = temp.*transpose(exp(pi*1i*(0:h.n-1)));
            temp = fftshift(temp);
            
            if h.M > h.N
                const = 1;
            else 
                const = 1/h.N;
            end

            % Multiplication with diagonal matrix
            temp = const .* transpose(1/h.n*1./infft.phi_hat(h.window,h.m,h.N,h.n,-h.N/2:h.N/2-1)) .* temp(h.n/2+1-h.N/2:h.n/2+1+h.N/2-1);
            
            % Set approximations
            h.fcheck = temp;
            
            % Computation time
            h.times.t_trafo = toc;
            
            % Set flag
            h.trafo_done = true;
        end %function
        
        function trafo_toep(h)
        % Computation of an iNFFT via Gohberg-Semencul formula.
            % Perform an adjoint NFFT
            tic
            plan = nfft(1,h.N,h.M);
            plan.x = -h.y(:);
            plan.f = h.f;
            nfft_adjoint(plan);
            temp = plan.fhat;

            % Multiplication with inverse Toeplitz matrix
            h.fcheck = infft.ToeplitzInvMul(h.D.x0,h.D.D1,h.D.D2,h.D.D3,h.D.D4,temp);
            
            % Computation time
            h.times.t_trafo = toc;
            
            % Set flag
            h.trafo_done = true;
        end %function
        
        function adjoint_quadratic(h)
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
            temp = infft.compute_coeff(h.M,plan,h.y_storage,h.x,alpha);
            
            % Final computation
            temp = h.d.*(temp + 1i*sum(alpha));
            
            % Set flag
            h.trafo_done = true;
            
            % If y got sorted, use the inverse permutation to get back the order of the user
            if h.perm
                h.ftilde(h.perm) = temp;
            else 
                h.ftilde = temp;
            end
            
            % Computation time
            h.times.t_trafo = toc;
        end %function
        
        function adjoint_rectangular(h)
        % Computation of an adjoint iNFFT for the rectangular setting.
            tic
            % Perform an NFFT
            % Multiplication with diagonal matrix
            temp = transpose(1./infft.phi_hat(h.window,h.m,h.N,h.n,-h.N/2:h.N/2-1)) .* h.fhat;
            
            % Perform an iFFT
            temp = ifft([zeros((h.n-h.N)/2,1);temp;zeros((h.n-h.N)/2,1)]);
            temp = temp.*transpose(exp(pi*1i*(0:h.n-1)));
            temp = fftshift(temp);
            
            % Multiplication with adjoint of optimized sparse matrix
            temp = h.B_opt'*temp;
            
            % Multiplication with specific constant
            if h.M > h.N
                const = 1;
            else 
                const = 1/h.N;
            end
            temp = const.*temp;

            % Set approximations
            h.ftilde = temp;

            % Computation time
            h.times.t_trafo = toc;
            
            % Set flag
            h.trafo_done = true;
        end %function

    end %methods
    
    
    % Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static) % Can also be used outside class by calling 'infft.methodsname(args)'
        
        function [vec,dist] = compute_sign(x,y,M)
        % Auxiliary function for computing the correction of sign for coefficients c.
        % Call from outside class: infft.compute_sign(x,y,M)
            count = zeros(M,1);
            i = 1;
            j = 1;

            if isempty(intersect(x,y))
                dist = 0; % Remember that no shift has to be done
            else
                dist = 1;
            end

            % Go through the nodes successively and count number of y's smaller than each node x
            while (i <= M && j <= M)
            if x(i) < y(j)
                diff = y(j)-x(i);
                if diff < dist
                    dist = diff; % Remember the smallest positive distance
                end
                i = i + 1;
                if i <= M
                    count(i) = count(i-1);
                end
            elseif x(i) > y(j)
                count(i) = count(i) + 1;
                j = j + 1;
            else % h.x(i) = h.y(j)
                if j < M
                    diff = y(j+1)-x(i);
                    if diff < dist
                        dist = diff;
                    end
                end
                i = i + 1;
                if i <= M
                    count(i) = count(i-1);
                end
                count(i-1) = count(i-1) + 1;
            end
            end

            vec = (-1).^count; % Vector with correct signs
        end%function
        
        function coeff = compute_coeff(M,plan,a,b,alpha)
        % Auxiliary function for computing needed coefficients.
        % Call from outside class: infft.compute_coeff(M,plan,a,b,alpha)
            coeff = zeros(M,1);

            % New index sets for nodes a
            A1 = a<-0.25;
            A2 = a>=0.25;
            A3 = and(not(A1),not(A2));

            % New index sets for nodes b
            B1 = b<0.25;
            B2 = not(B1);
            B3 = b>=-0.25;
            B4 = not(B3);
            %------------------------------------------------------------------
            % Case 1 (a < -0.25)
            % Computation without shift of nodes
            plan.x = b(B1)';
            plan.alpha = alpha(B1);
            plan.y = a(A1)';
            fastsum_trafo(plan)
            coeff(A1) = plan.f;

            % Computation with a shift
            plan.x = (b(B2)-1/4)';
            plan.alpha = alpha(B2);
            plan.y = (a(A1)+3/4)';
            fastsum_trafo(plan)
            coeff(A1) = coeff(A1) + plan.f;          
            %------------------------------------------------------------------
            % Case 2 (a >= 0.25)
            % Computation without shift of nodes
            plan.x = b(B3)';
            plan.alpha = alpha(B3);
            plan.y = a(A2)';
            fastsum_trafo(plan)
            coeff(A2) = plan.f;

            % Computation with a shift
            plan.x = (b(B4)+1/4)';
            plan.alpha = alpha(B4);
            plan.y = (a(A2)-3/4)';
            fastsum_trafo(plan)
            coeff(A2) = coeff(A2) + plan.f;
            %------------------------------------------------------------------
            % Case 3 (-0.25 <= a < 0.25)
            % Computation without shift of nodes
            plan.x = b';
            plan.alpha = alpha;
            plan.y = a(A3)';
            fastsum_trafo(plan)
            coeff(A3) = plan.f;
        end%function
        
        function phihat = phi_hat(window,m,N,n,v)
        % Computation of the Fourier transformed window function at points v
        % Call from outside class: infft.phi_hat(window,m,N,n,v)
            o = n/N;
            switch window
                case 'BSpline' % The Fourier transform of BSpline is the sinc function.
                    phihat = 1/n.*(infft.sinc(pi.*v/n)).^(2*m);
                case 'Gaussian'
                    b = 2*o/(2*o-1)*m/pi;
                    phihat = 1/n*exp(-b*(pi*v/n).^2);
                case 'Sinc' % The Fourier transform of the sinc function is BSpline.
                    phihat = infft.cardinal_bspline(2*m*v./((2*o-1)*N),2*m);
                case 'Kaiser'
                    phihat = zeros(1,length(v));
                    b = pi*(2-1/o);
                    ind = (abs(v./(1-1/(2*o))) <= n);
                    phihat(ind) = 1/n*besseli(0,m*sqrt(b^2-(2*pi*v(ind)./n).^2));
                case 'Bessel'
                    phihat = zeros(1,length(v));
                    b = pi*(2-1/o);
                    ind = abs(v) <= n*b/(2*pi);
                    phihat(ind) = 1/n*sinh(m*sqrt(b^2-4*pi^2*v(ind).^2/n^2))./sqrt(b^2-4*pi^2*v(ind).^2/n^2);
                    phihat(not(ind)) = 1/n*m*sinc(m*sqrt(4*pi^2*v(not(ind)).^2/n^2-b^2));
                case 'Dirichlet'
                    phihat = ones(size(v));
            end
        end%function
        
        function z = sinc(v)
        % Evaluation of the sinc function at points v
        % Call from outside class: infft.sinc(v)
            z = zeros(size(v));
            ind = (v == 0);
            z(ind) = 1;
            z(~ind) = sin(v(~ind))./v(~ind);
        end%function
        
        function z = cardinal_bspline(v,order)
        % Evaluation of the centered cardinal B-spline of given order at points v
        % Call from outside class: infft.cardinal_bspline(v,order)
        % Author: Franziska Nestler
            if ( order<=0 || order>20 || mod(order,1)~=0 )
                error('The given order has to be a natural number smaller than 20.')
            end
            
            z = zeros(size(v));

            for l = 1:length(v)
                u = v(l);

                if u<-order/2
                    j=0;
                elseif u>=order/2
                    j=0;
                else
                    j=ceil(u+order/2);
                end

                switch order
                case 1
                switch j
                     case 1; w = 1.0000000000000000e+00;
                     otherwise; w = 0;
                end 
                case 2
                switch j
                     case 1; w = 1.0000000000000000e+00 + u*( 1.0000000000000000e+00);
                     case 2; w = 1.0000000000000000e+00 + u*( -1.0000000000000000e+00);
                     otherwise; w = 0;
                end 
                case 3
                switch j
                     case 1; w = 1.1250000000000000e+00 + u*( 1.5000000000000000e+00 + u*( 5.0000000000000000e-01));
                     case 2; w = 7.5000000000000000e-01 + u*( 0.0000000000000000e+00 + u*( -1.0000000000000000e+00));
                     case 3; w = 1.1250000000000000e+00 + u*( -1.5000000000000000e+00 + u*( 5.0000000000000000e-01));
                     otherwise; w = 0;
                end 
                case 4
                switch j
                     case 1; w = 1.3333333333333333e+00 + u*( 2.0000000000000000e+00 + u*( 1.0000000000000000e+00 + u*( 1.6666666666666666e-01)));
                     case 2; w = 6.6666666666666663e-01 + u*( 0.0000000000000000e+00 + u*( -1.0000000000000000e+00 + u*( -5.0000000000000000e-01)));
                     case 3; w = 6.6666666666666663e-01 + u*( 0.0000000000000000e+00 + u*( -1.0000000000000000e+00 + u*( 5.0000000000000000e-01)));
                     case 4; w = 1.3333333333333333e+00 + u*( -2.0000000000000000e+00 + u*( 1.0000000000000000e+00 + u*( -1.6666666666666666e-01)));
                     otherwise; w = 0;
                end 
                case 5
                switch j
                     case 1; w = 1.6276041666666665e+00 + u*( 2.6041666666666665e+00 + u*( 1.5625000000000000e+00 + u*( 4.1666666666666663e-01 + u*( 4.1666666666666664e-02))));
                     case 2; w = 5.7291666666666663e-01 + u*( -2.0833333333333343e-01 + u*( -1.2500000000000000e+00 + u*( -8.3333333333333337e-01 + u*( -1.6666666666666666e-01))));
                     case 3; w = 5.9895833333333326e-01 + u*( 0.0000000000000000e+00 + u*( -6.2500000000000000e-01 + u*( 0.0000000000000000e+00 + u*( 2.5000000000000000e-01))));
                     case 4; w = 5.7291666666666663e-01 + u*( 2.0833333333333343e-01 + u*( -1.2500000000000000e+00 + u*( 8.3333333333333337e-01 + u*( -1.6666666666666666e-01))));
                     case 5; w = 1.6276041666666665e+00 + u*( -2.6041666666666665e+00 + u*( 1.5625000000000000e+00 + u*( -4.1666666666666663e-01 + u*( 4.1666666666666664e-02))));
                     otherwise; w = 0;
                end 
                case 6
                switch j
                     case 1; w = 2.0249999999999999e+00 + u*( 3.3750000000000000e+00 + u*( 2.2500000000000000e+00 + u*( 7.5000000000000000e-01 + u*( 1.2500000000000000e-01 + u*( 8.3333333333333332e-03)))));
                     case 2; w = 4.2499999999999999e-01 + u*( -6.2500000000000033e-01 + u*( -1.7500000000000000e+00 + u*( -1.2500000000000000e+00 + u*( -3.7499999999999994e-01 + u*( -4.1666666666666664e-02)))));
                     case 3; w = 5.5000000000000004e-01 + u*( -1.1102230246251565e-16 + u*( -4.9999999999999983e-01 + u*( 0.0000000000000000e+00 + u*( 2.5000000000000006e-01 + u*( 8.3333333333333329e-02)))));
                     case 4; w = 5.5000000000000004e-01 + u*( 1.1102230246251565e-16 + u*( -5.0000000000000000e-01 + u*( 0.0000000000000000e+00 + u*( 2.5000000000000000e-01 + u*( -8.3333333333333329e-02)))));
                     case 5; w = 4.2499999999999999e-01 + u*( 6.2500000000000022e-01 + u*( -1.7500000000000000e+00 + u*( 1.2500000000000000e+00 + u*( -3.7500000000000000e-01 + u*( 4.1666666666666664e-02)))));
                     case 6; w = 2.0249999999999999e+00 + u*( -3.3750000000000000e+00 + u*( 2.2500000000000000e+00 + u*( -7.5000000000000000e-01 + u*( 1.2500000000000000e-01 + u*( -8.3333333333333332e-03)))));
                     otherwise; w = 0;
                end 
                case 7
                switch j
                     case 1; w = 2.5531467013888887e+00 + u*( 4.3768229166666659e+00 + u*( 3.1263020833333335e+00 + u*( 1.1909722222222221e+00 + u*( 2.5520833333333331e-01 + u*( 2.9166666666666664e-02 + u*( 1.3888888888888889e-03))))));
                     case 2; w = 1.7955729166666637e-01 + u*( -1.3197916666666669e+00 + u*( -2.5703125000000000e+00 + u*( -1.8472222222222221e+00 + u*( -6.5625000000000000e-01 + u*( -1.1666666666666665e-01 + u*( -8.3333333333333332e-03))))));
                     case 3; w = 5.1178385416666694e-01 + u*( 9.1145833333331847e-03 + u*( -3.5546874999999956e-01 + u*( 1.2152777777777778e-01 + u*( 3.2812500000000006e-01 + u*( 1.4583333333333334e-01 + u*( 2.0833333333333332e-02))))));
                     case 4; w = 5.1102430555555556e-01 + u*( -1.2952601953960160e-16 + u*( -4.0104166666666652e-01 + u*( -7.4014868308343765e-17 + u*( 1.4583333333333334e-01 + u*( 0.0000000000000000e+00 + u*( -2.7777777777777776e-02))))));
                     case 5; w = 5.1178385416666694e-01 + u*( -9.1145833333330373e-03 + u*( -3.5546875000000000e-01 + u*( -1.2152777777777778e-01 + u*( 3.2812500000000000e-01 + u*( -1.4583333333333334e-01 + u*( 2.0833333333333332e-02))))));
                     case 6; w = 1.7955729166666645e-01 + u*( 1.3197916666666669e+00 + u*( -2.5703125000000000e+00 + u*( 1.8472222222222221e+00 + u*( -6.5625000000000000e-01 + u*( 1.1666666666666665e-01 + u*( -8.3333333333333332e-03))))));
                     case 7; w = 2.5531467013888887e+00 + u*( -4.3768229166666659e+00 + u*( 3.1263020833333335e+00 + u*( -1.1909722222222221e+00 + u*( 2.5520833333333331e-01 + u*( -2.9166666666666664e-02 + u*( 1.3888888888888889e-03))))));
                     otherwise; w = 0;
                end 
                case 8
                switch j
                     case 1; w = 3.2507936507936499e+00 + u*( 5.6888888888888891e+00 + u*( 4.2666666666666666e+00 + u*( 1.7777777777777775e+00 + u*( 4.4444444444444448e-01 + u*( 6.6666666666666666e-02 + u*( 5.5555555555555558e-03 + u*( 1.9841269841269841e-04)))))));
                     case 2; w = -2.2063492063492082e-01 + u*( -2.4111111111111110e+00 + u*( -3.8333333333333330e+00 + u*( -2.7222222222222219e+00 + u*( -1.0555555555555556e+00 + u*( -2.3333333333333334e-01 + u*( -2.7777777777777780e-02 + u*( -1.3888888888888889e-03)))))));
                     case 3; w = 4.9047619047619057e-01 + u*( 7.7777777777777932e-02 + u*( -9.9999999999999520e-02 + u*( 3.8888888888888917e-01 + u*( 5.0000000000000000e-01 + u*( 2.3333333333333336e-01 + u*( 4.9999999999999996e-02 + u*( 4.1666666666666666e-03)))))));
                     case 4; w = 4.7936507936507977e-01 + u*( -4.1930744590753681e-16 + u*( -3.3333333333333293e-01 + u*( -1.3481279584734043e-16 + u*( 1.1111111111111112e-01 + u*( 0.0000000000000000e+00 + u*( -2.7777777777777783e-02 + u*( -6.9444444444444432e-03)))))));
                     case 5; w = 4.7936507936507972e-01 + u*( -9.9127055770103257e-19 + u*( -3.3333333333333298e-01 + u*( -8.7231809077690869e-17 + u*( 1.1111111111111109e-01 + u*( -1.5860328923216521e-17 + u*( -2.7777777777777780e-02 + u*( 6.9444444444444432e-03)))))));
                     case 6; w = 4.9047619047619057e-01 + u*( -7.7777777777777682e-02 + u*( -1.0000000000000012e-01 + u*( -3.8888888888888901e-01 + u*( 5.0000000000000000e-01 + u*( -2.3333333333333336e-01 + u*( 4.9999999999999996e-02 + u*( -4.1666666666666666e-03)))))));
                     case 7; w = -2.2063492063492077e-01 + u*( 2.4111111111111114e+00 + u*( -3.8333333333333330e+00 + u*( 2.7222222222222219e+00 + u*( -1.0555555555555556e+00 + u*( 2.3333333333333334e-01 + u*( -2.7777777777777773e-02 + u*( 1.3888888888888889e-03)))))));
                     case 8; w = 3.2507936507936499e+00 + u*( -5.6888888888888891e+00 + u*( 4.2666666666666666e+00 + u*( -1.7777777777777775e+00 + u*( 4.4444444444444448e-01 + u*( -6.6666666666666666e-02 + u*( 5.5555555555555558e-03 + u*( -1.9841269841269841e-04)))))));
                     otherwise; w = 0;
                end 
                case 9
                switch j
                     case 1; w = 4.1704180036272316e+00 + u*( 7.4140764508928569e+00 + u*( 5.7665039062499996e+00 + u*( 2.5628906249999996e+00 + u*( 7.1191406250000011e-01 + u*( 1.2656250000000002e-01 + u*( 1.4062500000000000e-02 + u*( 8.9285714285714294e-04 + u*( 2.4801587301587302e-05))))))));
                     case 2; w = -8.5608956473214359e-01 + u*( -4.0750837053571427e+00 + u*( -5.7226562500000000e+00 + u*( -4.0023437500000014e+00 + u*( -1.6328125000000000e+00 + u*( -4.0937499999999993e-01 + u*( -6.2499999999999993e-02 + u*( -5.3571428571428572e-03 + u*( -1.9841269841269841e-04))))))));
                     case 3; w = 5.0630231584821439e-01 + u*( 2.8457031250000059e-01 + u*( 3.8085937500000017e-01 + u*( 8.8046875000000069e-01 + u*( 8.0859374999999989e-01 + u*( 3.7187500000000001e-01 + u*( 9.3749999999999986e-02 + u*( 1.2499999999999999e-02 + u*( 6.9444444444444447e-04))))))));
                     case 4; w = 4.5290876116071449e-01 + u*( -1.9531250000014036e-04 + u*( -2.8359374999999992e-01 + u*( -5.4687499999998678e-03 + u*( 7.0312499999999861e-02 + u*( -2.1874999999999967e-02 + u*( -3.7500000000000019e-02 + u*( -1.2500000000000001e-02 + u*( -1.3888888888888887e-03))))))));
                     case 5; w = 4.5292096819196492e-01 + u*( -4.1199682554449168e-16 + u*( -2.8222656249999944e-01 + u*( -2.1510571102112408e-16 + u*( 8.3984374999999944e-02 + u*( -1.3877787807814457e-17 + u*( -1.5625000000000014e-02 + u*( 8.6736173798840355e-19 + u*( 1.7361111111111108e-03))))))));
                     case 6; w = 4.5290876116071466e-01 + u*( 1.9531249999965074e-04 + u*( -2.8359374999999992e-01 + u*( 5.4687499999997776e-03 + u*( 7.0312499999999875e-02 + u*( 2.1874999999999933e-02 + u*( -3.7500000000000006e-02 + u*( 1.2499999999999999e-02 + u*( -1.3888888888888887e-03))))))));
                     case 7; w = 5.0630231584821439e-01 + u*( -2.8457031249999992e-01 + u*( 3.8085937500000017e-01 + u*( -8.8046875000000024e-01 + u*( 8.0859374999999989e-01 + u*( -3.7187500000000001e-01 + u*( 9.3750000000000000e-02 + u*( -1.2499999999999997e-02 + u*( 6.9444444444444447e-04))))))));
                     case 8; w = -8.5608956473214382e-01 + u*( 4.0750837053571427e+00 + u*( -5.7226562500000000e+00 + u*( 4.0023437500000005e+00 + u*( -1.6328125000000000e+00 + u*( 4.0937499999999993e-01 + u*( -6.2499999999999986e-02 + u*( 5.3571428571428563e-03 + u*( -1.9841269841269841e-04))))))));
                     case 9; w = 4.1704180036272316e+00 + u*( -7.4140764508928569e+00 + u*( 5.7665039062499996e+00 + u*( -2.5628906249999996e+00 + u*( 7.1191406250000011e-01 + u*( -1.2656250000000002e-01 + u*( 1.4062500000000000e-02 + u*( -8.9285714285714294e-04 + u*( 2.4801587301587302e-05))))))));
                     otherwise; w = 0;
                end 
                case 10
                switch j
                     case 1; w = 5.3822889109347445e+00 + u*( 9.6881200396825395e+00 + u*( 7.7504960317460307e+00 + u*( 3.6168981481481475e+00 + u*( 1.0850694444444446e+00 + u*( 2.1701388888888895e-01 + u*( 2.8935185185185189e-02 + u*( 2.4801587301587300e-03 + u*( 1.2400793650793653e-04 + u*( 2.7557319223985893e-06)))))))));
                     case 2; w = -1.8416969797178138e+00 + u*( -6.5658482142857135e+00 + u*( -8.5034722222222232e+00 + u*( -5.8645833333333339e+00 + u*( -2.4704861111111116e+00 + u*( -6.7187500000000000e-01 + u*( -1.1921296296296295e-01 + u*( -1.3392857142857140e-02 + u*( -8.6805555555555572e-04 + u*( -2.4801587301587302e-05)))))))));
                     case 3; w = 5.9915123456790131e-01 + u*( 7.5669642857142883e-01 + u*( 1.2599206349206369e+00 + u*( 1.7291666666666661e+00 + u*( 1.3263888888888888e+00 + u*( 5.9375000000000000e-01 + u*( 1.6203703703703703e-01 + u*( 2.6785714285714281e-02 + u*( 2.4801587301587300e-03 + u*( 9.9206349206349220e-05)))))))));
                     case 4; w = 4.2983906525573179e-01 + u*( -5.2083333333328057e-03 + u*( -2.6388888888888967e-01 + u*( -4.8611111111110529e-02 + u*( -6.9444444444446956e-03 + u*( -7.2916666666666644e-02 + u*( -6.0185185185185203e-02 + u*( -2.0833333333333336e-02 + u*( -3.4722222222222220e-03 + u*( -2.3148148148148149e-04)))))))));
                     case 5; w = 4.3041776895943618e-01 + u*( -2.4454782334950824e-18 + u*( -2.4305555555555569e-01 + u*( 1.3010426069826053e-16 + u*( 6.5972222222221863e-02 + u*( 2.8526563827174158e-17 + u*( -1.1574074074074106e-02 + u*( -1.5419764230904951e-18 + u*( 1.7361111111111110e-03 + u*( 3.4722222222222213e-04)))))))));
                     case 6; w = 4.3041776895943623e-01 + u*( -9.9439409253128842e-16 + u*( -2.4305555555555530e-01 + u*( -5.7168775886080107e-16 + u*( 6.5972222222221932e-02 + u*( -5.7053127654348317e-17 + u*( -1.1574074074074105e-02 + u*( 3.0839528461809902e-18 + u*( 1.7361111111111106e-03 + u*( -3.4722222222222213e-04)))))))));
                     case 7; w = 4.2983906525573234e-01 + u*( 5.2083333333329627e-03 + u*( -2.6388888888888878e-01 + u*( 4.8611111111110640e-02 + u*( -6.9444444444445291e-03 + u*( 7.2916666666666644e-02 + u*( -6.0185185185185182e-02 + u*( 2.0833333333333336e-02 + u*( -3.4722222222222216e-03 + u*( 2.3148148148148149e-04)))))))));
                     case 8; w = 5.9915123456790087e-01 + u*( -7.5669642857142860e-01 + u*( 1.2599206349206351e+00 + u*( -1.7291666666666665e+00 + u*( 1.3263888888888888e+00 + u*( -5.9375000000000000e-01 + u*( 1.6203703703703703e-01 + u*( -2.6785714285714288e-02 + u*( 2.4801587301587300e-03 + u*( -9.9206349206349220e-05)))))))));
                     case 9; w = -1.8416969797178138e+00 + u*( 6.5658482142857135e+00 + u*( -8.5034722222222232e+00 + u*( 5.8645833333333330e+00 + u*( -2.4704861111111112e+00 + u*( 6.7187499999999989e-01 + u*( -1.1921296296296295e-01 + u*( 1.3392857142857140e-02 + u*( -8.6805555555555551e-04 + u*( 2.4801587301587302e-05)))))))));
                     case 10; w = 5.3822889109347445e+00 + u*( -9.6881200396825395e+00 + u*( 7.7504960317460307e+00 + u*( -3.6168981481481475e+00 + u*( 1.0850694444444446e+00 + u*( -2.1701388888888895e-01 + u*( 2.8935185185185189e-02 + u*( -2.4801587301587300e-03 + u*( 1.2400793650793653e-04 + u*( -2.7557319223985893e-06)))))))));
                     otherwise; w = 0;
                end 
                case 11
                switch j
                     case 1; w = 6.9801356404084149e+00 + u*( 1.2691155709833485e+01 + u*( 1.0383672853500119e+01 + u*( 5.0345080501818780e+00 + u*( 1.6018889250578705e+00 + u*( 3.4950303819444450e-01 + u*( 5.2955005787037046e-02 + u*( 5.5018187830687838e-03 + u*( 3.7512400793650796e-04 + u*( 1.5156525573192242e-05 + u*( 2.7557319223985894e-07))))))))));
                     case 2; w = -3.3416489185689819e+00 + u*( -1.0246143310116292e+01 + u*( -1.2553626166449654e+01 + u*( -8.5579654431216952e+00 + u*( -3.6840729890046298e+00 + u*( -1.0600868055555557e+00 + u*( -2.0808015046296294e-01 + u*( -2.7645502645502647e-02 + u*( -2.3871527777777780e-03 + u*( -1.2125220458553793e-04 + u*( -2.7557319223985888e-06))))))))));
                     case 3; w = 8.3935043227315109e-01 + u*( 1.6995691208612389e+00 + u*( 2.8051469590928821e+00 + u*( 3.1439569382440475e+00 + u*( 2.1668882016782396e+00 + u*( 9.4595703124999986e-01 + u*( 2.6954933449074076e-01 + u*( 5.0334821428571430e-02 + u*( 5.9678819444444441e-03 + u*( 4.0922619047619049e-04 + u*( 1.2400793650793652e-05))))))))));
                     case 4; w = 4.0571875731991314e-01 + u*( -3.4957578951720865e-02 + u*( -3.1700110057043507e-01 + u*( -1.8633432539682604e-01 + u*( -1.6431568287037052e-01 + u*( -1.7302083333333337e-01 + u*( -1.0344328703703705e-01 + u*( -3.4920634920634935e-02 + u*( -6.8204365079365071e-03 + u*( -7.2751322751322749e-04 + u*( -3.3068783068783071e-05))))))))));
                     case 5; w = 4.1096276715529234e-01 + u*( 2.4866174778753639e-06 + u*( -2.1212090386284888e-01 + u*( 1.1935763888936505e-04 + u*( 5.3213614004629026e-02 + u*( 1.0026041666666599e-03 + u*( -6.7635995370370679e-03 + u*( 1.9097222222222133e-03 + u*( 2.3871527777777780e-03 + u*( 6.3657407407407391e-04 + u*( 5.7870370370370366e-05))))))))));
                     case 6; w = 4.1096264282441891e-01 + u*( -5.5155532918682580e-16 + u*( -2.1214328342014022e-01 + u*( -2.1510571102112408e-16 + u*( 5.2795862268517693e-02 + u*( -6.9388939039072284e-18 + u*( -8.4346064814815325e-03 + u*( 6.9388939039072288e-19 + u*( 9.5486111111111119e-04 + u*( -4.3368086899420180e-20 + u*( -6.9444444444444431e-05))))))))));
                     case 7; w = 4.1096276715529312e-01 + u*( -2.4866174781681637e-06 + u*( -2.1212090386284718e-01 + u*( -1.1935763888944589e-04 + u*( 5.3213614004629338e-02 + u*( -1.0026041666666237e-03 + u*( -6.7635995370370653e-03 + u*( -1.9097222222222091e-03 + u*( 2.3871527777777771e-03 + u*( -6.3657407407407391e-04 + u*( 5.7870370370370366e-05))))))))));
                     case 8; w = 4.0571875731991297e-01 + u*( 3.4957578951719165e-02 + u*( -3.1700110057043585e-01 + u*( 1.8633432539682518e-01 + u*( -1.6431568287037035e-01 + u*( 1.7302083333333326e-01 + u*( -1.0344328703703705e-01 + u*( 3.4920634920634921e-02 + u*( -6.8204365079365089e-03 + u*( 7.2751322751322749e-04 + u*( -3.3068783068783071e-05))))))))));
                     case 9; w = 8.3935043227315076e-01 + u*( -1.6995691208612371e+00 + u*( 2.8051469590928817e+00 + u*( -3.1439569382440480e+00 + u*( 2.1668882016782400e+00 + u*( -9.4595703124999986e-01 + u*( 2.6954933449074076e-01 + u*( -5.0334821428571430e-02 + u*( 5.9678819444444441e-03 + u*( -4.0922619047619049e-04 + u*( 1.2400793650793652e-05))))))))));
                     case 10; w = -3.3416489185689819e+00 + u*( 1.0246143310116292e+01 + u*( -1.2553626166449654e+01 + u*( 8.5579654431216934e+00 + u*( -3.6840729890046298e+00 + u*( 1.0600868055555555e+00 + u*( -2.0808015046296294e-01 + u*( 2.7645502645502647e-02 + u*( -2.3871527777777780e-03 + u*( 1.2125220458553791e-04 + u*( -2.7557319223985888e-06))))))))));
                     case 11; w = 6.9801356404084149e+00 + u*( -1.2691155709833485e+01 + u*( 1.0383672853500119e+01 + u*( -5.0345080501818780e+00 + u*( 1.6018889250578705e+00 + u*( -3.4950303819444450e-01 + u*( 5.2955005787037046e-02 + u*( -5.5018187830687838e-03 + u*( 3.7512400793650796e-04 + u*( -1.5156525573192242e-05 + u*( 2.7557319223985894e-07))))))))));
                     otherwise; w = 0;
                end 
                case 12
                switch j
                     case 1; w = 9.0888311688311667e+00 + u*( 1.6662857142857138e+01 + u*( 1.3885714285714283e+01 + u*( 6.9428571428571413e+00 + u*( 2.3142857142857149e+00 + u*( 5.4000000000000015e-01 + u*( 9.0000000000000024e-02 + u*( 1.0714285714285716e-02 + u*( 8.9285714285714305e-04 + u*( 4.9603174603174603e-05 + u*( 1.6534391534391537e-06 + u*( 2.5052108385441720e-08)))))))))));
                     case 2; w = -5.5901385882635886e+00 + u*( -1.5630876322751321e+01 + u*( -1.8408019179894186e+01 + u*( -1.2433382936507940e+01 + u*( -5.4362103174603176e+00 + u*( -1.6301388888888892e+00 + u*( -3.4402777777777782e-01 + u*( -5.1289682539682538e-02 + u*( -5.3075396825396828e-03 + u*( -3.6375661375661386e-04 + u*( -1.4880952380952381e-05 + u*( -2.7557319223985888e-07)))))))))));
                     case 3; w = 1.3448878667628679e+00 + u*( 3.4404464285714282e+00 + u*( 5.4311342592592586e+00 + u*( 5.4459821428571402e+00 + u*( 3.5034722222222210e+00 + u*( 1.4987499999999996e+00 + u*( 4.3819444444444450e-01 + u*( 8.8392857142857148e-02 + u*( 1.2152777777777778e-02 + u*( 1.0912698412698413e-03 + u*( 5.7870370370370373e-05 + u*( 1.3778659611992946e-06)))))))))));
                     case 4; w = 3.6854858104857957e-01 + u*( -1.3946428571428326e-01 + u*( -5.3538359788359824e-01 + u*( -5.2053571428571310e-01 + u*( -4.7420634920634980e-01 + u*( -3.5750000000000021e-01 + u*( -1.8055555555555561e-01 + u*( -5.8928571428571448e-02 + u*( -1.2400793650793652e-02 + u*( -1.6369047619047622e-03 + u*( -1.2400793650793650e-04 + u*( -4.1335978835978839e-06)))))))))));
                     case 5; w = 3.9394540644540721e-01 + u*( 2.1825396825182083e-04 + u*( -1.8617724867724822e-01 + u*( 3.2738095238083755e-03 + u*( 4.9603174603174205e-02 + u*( 9.1666666666664915e-03 + u*( 2.7777777777777675e-03 + u*( 6.5476190476190192e-03 + u*( 3.9682539682539689e-03 + u*( 1.0912698412698409e-03 + u*( 1.4880952380952380e-04 + u*( 8.2671957671957678e-06)))))))))));
                     case 6; w = 3.9392556517556448e-01 + u*( 5.6286763900459553e-16 + u*( -1.8726851851852244e-01 + u*( -1.3273098687069132e-16 + u*( 4.3055555555554292e-02 + u*( -9.8997514731312790e-17 + u*( -6.3888888888889292e-03 + u*( -1.1039149392579681e-17 + u*( 6.9444444444444621e-04 + u*( -1.9712766772463716e-19 + u*( -6.9444444444444431e-05 + u*( -1.1574074074074073e-05)))))))))));
                     case 7; w = 3.9392556517556537e-01 + u*( -4.6512889223589778e-16 + u*( -1.8726851851852064e-01 + u*( 3.6751989551412045e-16 + u*( 4.3055555555554653e-02 + u*( 1.6531126215388073e-16 + u*( -6.3888888888889214e-03 + u*( 1.2616170734376778e-17 + u*( 6.9444444444444556e-04 + u*( -3.9425533544927435e-20 + u*( -6.9444444444444417e-05 + u*( 1.1574074074074073e-05)))))))))));
                     case 8; w = 3.9394540644540721e-01 + u*( -2.1825396825495211e-04 + u*( -1.8617724867724822e-01 + u*( -3.2738095238095664e-03 + u*( 4.9603174603174142e-02 + u*( -9.1666666666665921e-03 + u*( 2.7777777777776894e-03 + u*( -6.5476190476190374e-03 + u*( 3.9682539682539646e-03 + u*( -1.0912698412698411e-03 + u*( 1.4880952380952380e-04 + u*( -8.2671957671957678e-06)))))))))));
                     case 9; w = 3.6854858104858013e-01 + u*( 1.3946428571428304e-01 + u*( -5.3538359788359802e-01 + u*( 5.2053571428571266e-01 + u*( -4.7420634920634908e-01 + u*( 3.5750000000000010e-01 + u*( -1.8055555555555555e-01 + u*( 5.8928571428571434e-02 + u*( -1.2400793650793652e-02 + u*( 1.6369047619047617e-03 + u*( -1.2400793650793650e-04 + u*( 4.1335978835978839e-06)))))))))));
                     case 10; w = 1.3448878667628670e+00 + u*( -3.4404464285714287e+00 + u*( 5.4311342592592586e+00 + u*( -5.4459821428571420e+00 + u*( 3.5034722222222214e+00 + u*( -1.4987499999999998e+00 + u*( 4.3819444444444461e-01 + u*( -8.8392857142857148e-02 + u*( 1.2152777777777775e-02 + u*( -1.0912698412698413e-03 + u*( 5.7870370370370373e-05 + u*( -1.3778659611992946e-06)))))))))));
                     case 11; w = -5.5901385882635886e+00 + u*( 1.5630876322751321e+01 + u*( -1.8408019179894186e+01 + u*( 1.2433382936507938e+01 + u*( -5.4362103174603176e+00 + u*( 1.6301388888888888e+00 + u*( -3.4402777777777782e-01 + u*( 5.1289682539682538e-02 + u*( -5.3075396825396819e-03 + u*( 3.6375661375661375e-04 + u*( -1.4880952380952383e-05 + u*( 2.7557319223985888e-07)))))))))));
                     case 12; w = 9.0888311688311667e+00 + u*( -1.6662857142857138e+01 + u*( 1.3885714285714283e+01 + u*( -6.9428571428571413e+00 + u*( 2.3142857142857149e+00 + u*( -5.4000000000000015e-01 + u*( 9.0000000000000024e-02 + u*( -1.0714285714285716e-02 + u*( 8.9285714285714305e-04 + u*( -4.9603174603174603e-05 + u*( 1.6534391534391537e-06 + u*( -2.5052108385441720e-08)))))))))));
                     otherwise; w = 0;
                end 
                case 13
                switch j
                     case 1; w = 1.1874718295524923e+01 + u*( 2.1922556853276788e+01 + u*( 1.8549855798926512e+01 + u*( 9.5127465635520565e+00 + u*( 3.2928738104603297e+00 + u*( 8.1055355334408086e-01 + u*( 1.4548397111304015e-01 + u*( 1.9184699487433867e-02 + u*( 1.8446826430224871e-03 + u*( 1.2613214653145211e-04 + u*( 5.8214836860670193e-06 + u*( 1.6283870450537117e-07 + u*( 2.0876756987868100e-09))))))))))));
                     case 2; w = -8.9202691331918125e+00 + u*( -2.3448324809377926e+01 + u*( -2.6821025863728199e+01 + u*( -1.7984757474420494e+01 + u*( -7.9561051141648065e+00 + u*( -2.4618766792741402e+00 + u*( -5.4866789641203706e-01 + u*( -8.8994812334656115e-02 + u*( -1.0448443700396826e-02 + u*( -8.6725180041152285e-04 + u*( -4.8363095238095247e-05 + u*( -1.6283870450537117e-06 + u*( -2.5052108385441717e-08))))))))))));
                     case 3; w = 2.3078704748693437e+00 + u*( 6.4933808121185033e+00 + u*( 9.7743921181007654e+00 + u*( 9.1229595491565103e+00 + u*( 5.5977533976236957e+00 + u*( 2.3572730138062163e+00 + u*( 7.0074128327546292e-01 + u*( 1.4898788855820111e-01 + u*( 2.2604709201388889e-02 + u*( 2.3972571281599057e-03 + u*( 1.6927083333333337e-04 + u*( 7.1649029982363331e-06 + u*( 1.3778659611992947e-07))))))))))));
                     case 4; w = 2.9022154571295156e-01 + u*( -4.2427265927484553e-01 + u*( -1.0962061940887806e+00 + u*( -1.2299912243573414e+00 + u*( -1.0577149567780679e+00 + u*( -6.8522680534887614e-01 + u*( -3.1342532310956794e-01 + u*( -9.9379443617724925e-02 + u*( -2.1746600115740752e-02 + u*( -3.2346551660787782e-03 + u*( -3.1346450617283954e-04 + u*( -1.7912257495590829e-05 + u*( -4.5928865373309819e-07))))))))))));
                     case 5; w = 3.7919269114211729e-01 + u*( 2.7888387851622309e-03 + u*( -1.5667089835676884e-01 + u*( 2.2722503285348961e-02 + u*( 6.9727398100352347e-02 + u*( 3.6336301773313151e-02 + u*( 2.3304126880786955e-02 + u*( 1.6070653521825361e-02 + u*( 7.1159241691468233e-03 + u*( 1.8964602623456788e-03 + u*( 3.0226934523809519e-04 + u*( 2.6868386243386244e-05 + u*( 1.0333994708994710e-06))))))))))));
                     case 6; w = 3.7884408367010836e-01 + u*( -2.0990931083435783e-08 + u*( -1.6689671753575808e-01 + u*( -1.5393346312846418e-06 + u*( 3.5641334170385862e-02 + u*( -2.2166418651253612e-05 + u*( -4.9746817129629447e-03 + u*( -8.8665674603212124e-05 + u*( 3.8287450396825841e-04 + u*( -9.8517416225749903e-05 + u*( -9.6726190476190466e-05 + u*( -2.1494708994708992e-05 + u*( -1.6534391534391535e-06))))))))))));
                     case 7; w = 3.7884408454472845e-01 + u*( 1.3561879811528072e-15 + u*( -1.6689648663556647e-01 + u*( 8.4999643319242735e-16 + u*( 3.5648261176214106e-02 + u*( 9.6710833785706996e-17 + u*( -4.9229600694444126e-03 + u*( -9.3964188282077045e-19 + u*( 4.9370659722222730e-04 + u*( -2.1684043449710089e-19 + u*( -3.7615740740740696e-05 + u*( 2.2587545260114674e-21 + u*( 1.9290123456790124e-06))))))))))));
                     case 8; w = 3.7884408367010813e-01 + u*( 2.0990927702110751e-08 + u*( -1.6689671753575838e-01 + u*( 1.5393346297883930e-06 + u*( 3.5641334170385473e-02 + u*( 2.2166418650956720e-05 + u*( -4.9746817129630817e-03 + u*( 8.8665674603176291e-05 + u*( 3.8287450396825234e-04 + u*( 9.8517416225749226e-05 + u*( -9.6726190476190466e-05 + u*( 2.1494708994708992e-05 + u*( -1.6534391534391535e-06))))))))))));
                     case 9; w = 3.7919269114211857e-01 + u*( -2.7888387851645394e-03 + u*( -1.5667089835676726e-01 + u*( -2.2722503285349218e-02 + u*( 6.9727398100352833e-02 + u*( -3.6336301773313096e-02 + u*( 2.3304126880786948e-02 + u*( -1.6070653521825368e-02 + u*( 7.1159241691468198e-03 + u*( -1.8964602623456788e-03 + u*( 3.0226934523809519e-04 + u*( -2.6868386243386247e-05 + u*( 1.0333994708994710e-06))))))))))));
                     case 10; w = 2.9022154571295128e-01 + u*( 4.2427265927484364e-01 + u*( -1.0962061940887808e+00 + u*( 1.2299912243573401e+00 + u*( -1.0577149567780673e+00 + u*( 6.8522680534887614e-01 + u*( -3.1342532310956783e-01 + u*( 9.9379443617724869e-02 + u*( -2.1746600115740745e-02 + u*( 3.2346551660787774e-03 + u*( -3.1346450617283949e-04 + u*( 1.7912257495590826e-05 + u*( -4.5928865373309819e-07))))))))))));
                     case 11; w = 2.3078704748693428e+00 + u*( -6.4933808121185059e+00 + u*( 9.7743921181007636e+00 + u*( -9.1229595491565139e+00 + u*( 5.5977533976236975e+00 + u*( -2.3572730138062172e+00 + u*( 7.0074128327546303e-01 + u*( -1.4898788855820105e-01 + u*( 2.2604709201388889e-02 + u*( -2.3972571281599057e-03 + u*( 1.6927083333333337e-04 + u*( -7.1649029982363331e-06 + u*( 1.3778659611992947e-07))))))))))));
                     case 12; w = -8.9202691331918125e+00 + u*( 2.3448324809377926e+01 + u*( -2.6821025863728199e+01 + u*( 1.7984757474420494e+01 + u*( -7.9561051141648056e+00 + u*( 2.4618766792741398e+00 + u*( -5.4866789641203706e-01 + u*( 8.8994812334656115e-02 + u*( -1.0448443700396824e-02 + u*( 8.6725180041152263e-04 + u*( -4.8363095238095240e-05 + u*( 1.6283870450537117e-06 + u*( -2.5052108385441717e-08))))))))))));
                     case 13; w = 1.1874718295524923e+01 + u*( -2.1922556853276788e+01 + u*( 1.8549855798926512e+01 + u*( -9.5127465635520565e+00 + u*( 3.2928738104603297e+00 + u*( -8.1055355334408086e-01 + u*( 1.4548397111304015e-01 + u*( -1.9184699487433867e-02 + u*( 1.8446826430224871e-03 + u*( -1.2613214653145211e-04 + u*( 5.8214836860670193e-06 + u*( -1.6283870450537117e-07 + u*( 2.0876756987868100e-09))))))))))));
                     otherwise; w = 0;
                end 
                case 14
                switch j
                     case 1; w = 1.5559448654322788e+01 + u*( 2.8896118929456602e+01 + u*( 2.4768101939534233e+01 + u*( 1.2973767682613174e+01 + u*( 4.6334884580761324e+00 + u*( 1.1914684606481487e+00 + u*( 2.2694637345679020e-01 + u*( 3.2420910493827164e-02 + u*( 3.4736689814814817e-03 + u*( 2.7568801440329221e-04 + u*( 1.5753600823045270e-05 + u*( 6.1377665544332219e-07 + u*( 1.4613729891507670e-08 + u*( 1.6059043836821616e-10)))))))))))));
                     case 2; w = -1.3804467429593295e+01 + u*( -3.4725699252361579e+01 + u*( -3.8853716242283944e+01 + u*( -2.5906232317386841e+01 + u*( -1.1566511541923870e+01 + u*( -3.6685315393518532e+00 + u*( -8.5305362654320993e-01 + u*( -1.4757908950617285e-01 + u*( -1.9026331018518521e-02 + u*( -1.8076453189300417e-03 + u*( -1.2313528806584363e-04 + u*( -5.6993546576879917e-06 + u*( -1.6075102880658435e-07 + u*( -2.0876756987868096e-09)))))))))))));
                     case 3; w = 4.0345583168760264e+00 + u*( 1.1655767688458655e+01 + u*( 1.6804044086700323e+01 + u*( 1.4909458590534967e+01 + u*( 8.8413339120370367e+00 + u*( 3.6782928240740733e+00 + u*( 1.1060995370370370e+00 + u*( 2.4425154320987655e-01 + u*( 3.9748263888888899e-02 + u*( 4.7228652263374494e-03 + u*( 3.9930555555555563e-04 + u*( 2.2797418630751967e-05 + u*( 7.8914141414141425e-07 + u*( 1.2526054192720860e-08)))))))))))));
                     case 4; w = 1.1171506958833423e-01 + u*( -1.0934728652263375e+00 + u*( -2.3198167438271580e+00 + u*( -2.6207471707818968e+00 + u*( -2.1150446887860097e+00 + u*( -1.2520775462962976e+00 + u*( -5.3735725308642013e-01 + u*( -1.6661265432098765e-01 + u*( -3.7288773148148172e-02 + u*( -5.9767232510288086e-03 + u*( -6.7065329218106992e-04 + u*( -5.0154320987654331e-05 + u*( -2.2505144032921812e-06 + u*( -4.5928865373309820e-08)))))))))))));
                     case 5; w = 3.6800413208833282e-01 + u*( 1.7113072273659553e-02 + u*( -9.8644868827161739e-02 + u*( 9.4018454218104050e-02 + u*( 1.4725999871399059e-01 + u*( 1.0530526620370263e-01 + u*( 6.5923996913580168e-02 + u*( 3.4481095679012294e-02 + u*( 1.2984664351851844e-02 + u*( 3.3331725823045263e-03 + u*( 5.7066615226337438e-04 + u*( 6.2692901234567892e-05 + u*( 4.0187757201646094e-06 + u*( 1.1482216343327455e-07)))))))))))));
                     case 6; w = 3.6537038723236331e-01 + u*( -6.2692901262742653e-06 + u*( -1.5000289351852644e-01 + u*( -1.3792438271897567e-04 + u*( 2.9564525462960980e-02 + u*( -6.2065972222282459e-04 + u*( -4.6932870370370947e-03 + u*( -8.2754629629633357e-04 + u*( -2.5607638888888734e-04 + u*( -3.4481095679012323e-04 + u*( -1.6493055555555559e-04 + u*( -3.7615740740740737e-05 + u*( -4.3402777777777778e-06 + u*( -2.0667989417989419e-07)))))))))))));
                     case 7; w = 3.6537086948545300e-01 + u*( -4.2249810617422838e-15 + u*( -1.4996527777778026e-01 + u*( -2.1620807227687110e-15 + u*( 2.9909336419752904e-02 + u*( -3.7071374282292823e-16 + u*( -3.8657407407405330e-03 + u*( -3.0858061832279744e-17 + u*( 3.6458333333334490e-04 + u*( -3.0024060161137045e-19 + u*( -2.7006172839506100e-05 + u*( 7.2975146224985872e-21 + u*( 1.9290123456790124e-06 + u*( 2.7557319223985888e-07)))))))))))));
                     case 8; w = 3.6537086948545089e-01 + u*( 3.5866003099898934e-15 + u*( -1.4996527777778454e-01 + u*( 1.4151423104552569e-15 + u*( 2.9909336419751239e-02 + u*( -1.4594195243326033e-16 + u*( -3.8657407407407807e-03 + u*( -3.2534405191276562e-17 + u*( 3.6458333333333633e-04 + u*( -1.1342422727540662e-18 + u*( -2.7006172839506137e-05 + u*( -3.1275062667851091e-21 + u*( 1.9290123456790120e-06 + u*( -2.7557319223985888e-07)))))))))))));
                     case 9; w = 3.6537038723236431e-01 + u*( 6.2692901236420062e-06 + u*( -1.5000289351852428e-01 + u*( 1.3792438271880398e-04 + u*( 2.9564525462961413e-02 + u*( 6.2065972222283196e-04 + u*( -4.6932870370371416e-03 + u*( 8.2754629629630994e-04 + u*( -2.5607638888889141e-04 + u*( 3.4481095679012258e-04 + u*( -1.6493055555555556e-04 + u*( 3.7615740740740737e-05 + u*( -4.3402777777777778e-06 + u*( 2.0667989417989419e-07)))))))))))));
                     case 10; w = 3.6800413208833366e-01 + u*( -1.7113072273663724e-02 + u*( -9.8644868827161045e-02 + u*( -9.4018454218104675e-02 + u*( 1.4725999871399117e-01 + u*( -1.0530526620370259e-01 + u*( 6.5923996913580155e-02 + u*( -3.4481095679012308e-02 + u*( 1.2984664351851843e-02 + u*( -3.3331725823045263e-03 + u*( 5.7066615226337438e-04 + u*( -6.2692901234567892e-05 + u*( 4.0187757201646094e-06 + u*( -1.1482216343327455e-07)))))))))))));
                     case 11; w = 1.1171506958833281e-01 + u*( 1.0934728652263326e+00 + u*( -2.3198167438271593e+00 + u*( 2.6207471707818946e+00 + u*( -2.1150446887860097e+00 + u*( 1.2520775462962972e+00 + u*( -5.3735725308641968e-01 + u*( 1.6661265432098765e-01 + u*( -3.7288773148148151e-02 + u*( 5.9767232510288069e-03 + u*( -6.7065329218106992e-04 + u*( 5.0154320987654318e-05 + u*( -2.2505144032921812e-06 + u*( 4.5928865373309820e-08)))))))))))));
                     case 12; w = 4.0345583168760273e+00 + u*( -1.1655767688458658e+01 + u*( 1.6804044086700330e+01 + u*( -1.4909458590534971e+01 + u*( 8.8413339120370384e+00 + u*( -3.6782928240740751e+00 + u*( 1.1060995370370370e+00 + u*( -2.4425154320987655e-01 + u*( 3.9748263888888899e-02 + u*( -4.7228652263374485e-03 + u*( 3.9930555555555563e-04 + u*( -2.2797418630751967e-05 + u*( 7.8914141414141425e-07 + u*( -1.2526054192720860e-08)))))))))))));
                     case 13; w = -1.3804467429593299e+01 + u*( 3.4725699252361579e+01 + u*( -3.8853716242283951e+01 + u*( 2.5906232317386834e+01 + u*( -1.1566511541923868e+01 + u*( 3.6685315393518527e+00 + u*( -8.5305362654320993e-01 + u*( 1.4757908950617285e-01 + u*( -1.9026331018518521e-02 + u*( 1.8076453189300413e-03 + u*( -1.2313528806584363e-04 + u*( 5.6993546576879917e-06 + u*( -1.6075102880658435e-07 + u*( 2.0876756987868096e-09)))))))))))));
                     case 14; w = 1.5559448654322788e+01 + u*( -2.8896118929456602e+01 + u*( 2.4768101939534233e+01 + u*( -1.2973767682613174e+01 + u*( 4.6334884580761324e+00 + u*( -1.1914684606481487e+00 + u*( 2.2694637345679020e-01 + u*( -3.2420910493827164e-02 + u*( 3.4736689814814817e-03 + u*( -2.7568801440329221e-04 + u*( 1.5753600823045270e-05 + u*( -6.1377665544332219e-07 + u*( 1.4613729891507670e-08 + u*( -1.6059043836821616e-10)))))))))))));
                     otherwise; w = 0;
                end 
                case 15
                switch j
                     case 1; w = 2.0438514873693794e+01 + u*( 3.8151894430895062e+01 + u*( 3.3064975173442399e+01 + u*( 1.7634653425835946e+01 + u*( 6.4660395894731808e+00 + u*( 1.7242772238595152e+00 + u*( 3.4485544477190305e-01 + u*( 5.2549401108099490e-02 + u*( 6.1307634626116070e-03 + u*( 5.4495675223214287e-04 + u*( 3.6330450148809523e-05 + u*( 1.7614763708513713e-06 + u*( 5.8715879028379031e-08 + u*( 1.2044282877616211e-09 + u*( 1.1470745597729726e-11))))))))))))));
                     case 2; w = -2.0910950619651938e+01 + u*( -5.0908492785541867e+01 + u*( -5.5995412042994573e+01 + u*( -3.7171738707356027e+01 + u*( -1.6721280159184968e+01 + u*( -5.4102826988045321e+00 + u*( -1.3015814604582610e+00 + u*( -2.3693401080050075e-01 + u*( -3.2838157371238429e-02 + u*( -3.4518556409832463e-03 + u*( -2.7111665702160499e-04 + u*( -1.5438361792528462e-05 + u*( -6.0281635802469145e-07 + u*( -1.4453139453139453e-08 + u*( -1.6059043836821613e-10))))))))))))));
                     case 3; w = 7.0053971489285942e+00 + u*( 2.0151301534481288e+01 + u*( 2.7984344880669187e+01 + u*( 2.3904448146217611e+01 + u*( 1.3816813267601860e+01 + u*( 5.6944785472997692e+00 + u*( 1.7269897884792751e+00 + u*( 3.9237949547223405e-01 + u*( 6.7279900444878471e-02 + u*( 8.6836665185460772e-03 + u*( 8.3211263020833352e-04 + u*( 5.7502417528459207e-05 + u*( 2.7126736111111116e-06 + u*( 7.8287838704505374e-08 + u*( 1.0438378493934050e-09))))))))))))));
                     case 4; w = -2.8209730822648066e-01 + u*( -2.5209034433345012e+00 + u*( -4.7643956428425094e+00 + u*( -5.2055434302372499e+00 + u*( -3.9726260291205540e+00 + u*( -2.2119389179101931e+00 + u*( -9.0848269992404551e-01 + u*( -2.7694685078892683e-01 + u*( -6.2866889105902796e-02 + u*( -1.0597339340828927e-02 + u*( -1.3102213541666666e-03 + u*( -1.1561548019881353e-04 + u*( -6.9049873737373758e-06 + u*( -2.5052108385441723e-07 + u*( -4.1753513975736200e-09))))))))))))));
                     case 5; w = 3.6605477347366966e-01 + u*( 7.1704883466104596e-02 + u*( 5.0448392644339819e-02 + u*( 2.9713546746201114e-01 + u*( 3.5090739050028985e-01 + u*( 2.5865160758743305e-01 + u*( 1.5034181100350780e-01 + u*( 6.8791765024152032e-02 + u*( 2.3567764847366880e-02 + u*( 5.8664042693176788e-03 + u*( 1.0417420187114198e-03 + u*( 1.2874435074955906e-04 + u*( 1.0549286265432096e-05 + u*( 5.1669973544973548e-07 + u*( 1.1482216343327455e-08))))))))))))));
                     case 6; w = 3.5322239672907541e-01 + u*( -1.5642630361125927e-04 + u*( -1.3639101275692073e-01 + u*( -1.8075811799948417e-03 + u*( 2.2070036994084653e-02 + u*( -4.4182752175320316e-03 + u*( -7.5001186794705450e-03 + u*( -3.3645456880669446e-03 + u*( -1.6869439019097238e-03 + u*( -8.6818473048941736e-04 + u*( -3.0517578125000011e-04 + u*( -6.7170965608465618e-05 + u*( -9.0422453703703682e-06 + u*( -6.8893298059964710e-07 + u*( -2.2964432686654910e-08))))))))))))));
                     case 7; w = 3.5323915670369133e-01 + u*( 1.2614486916819811e-10 + u*( -1.3571316489467647e-01 + u*( 1.3119327370210617e-08 + u*( 2.5383959876165763e-02 + u*( 2.8862524260785272e-07 + u*( -3.0815548366968531e-03 + u*( 1.9791445223848725e-06 + u*( 2.7686225043403783e-04 + u*( 4.6180038855831713e-06 + u*( -1.4241536458333318e-05 + u*( 3.3585482804232957e-06 + u*( 2.7126736111111116e-06 + u*( 5.1669973544973538e-07 + u*( 3.4446649029982367e-08))))))))))))));
                     case 8; w = 3.5323915669918948e-01 + u*( -3.7627802471253198e-15 + u*( -1.3571316653458879e-01 + u*( -3.9208493232655638e-15 + u*( 2.5383887719855336e-02 + u*( -1.1909351001031243e-15 + u*( -3.0824207124253993e-03 + u*( -9.8495895219692053e-17 + u*( 2.7339874751985114e-04 + u*( -1.9360753080098294e-18 + u*( -1.8859540343915360e-05 + u*( -9.6803765400491462e-22 + u*( 1.0333994708994701e-06 + u*( 0.0000000000000000e+00 + u*( -3.9367598891408411e-08))))))))))))));
                     case 9; w = 3.5323915670368994e-01 + u*( -1.2614282685044184e-10 + u*( -1.3571316489467822e-01 + u*( -1.3119325577071029e-08 + u*( 2.5383959876165180e-02 + u*( -2.8862524270901942e-07 + u*( -3.0815548366969958e-03 + u*( -1.9791445224379912e-06 + u*( 2.7686225043403240e-04 + u*( -4.6180038855843724e-06 + u*( -1.4241536458333278e-05 + u*( -3.3585482804232868e-06 + u*( 2.7126736111111104e-06 + u*( -5.1669973544973538e-07 + u*( 3.4446649029982367e-08))))))))))))));
                     case 10; w = 3.5322239672907640e-01 + u*( 1.5642630360625928e-04 + u*( -1.3639101275691862e-01 + u*( 1.8075811799945099e-03 + u*( 2.2070036994085235e-02 + u*( 4.4182752175320377e-03 + u*( -7.5001186794705780e-03 + u*( 3.3645456880669099e-03 + u*( -1.6869439019097259e-03 + u*( 8.6818473048941682e-04 + u*( -3.0517578125000000e-04 + u*( 6.7170965608465605e-05 + u*( -9.0422453703703699e-06 + u*( 6.8893298059964721e-07 + u*( -2.2964432686654910e-08))))))))))))));
                     case 11; w = 3.6605477347366888e-01 + u*( -7.1704883466112257e-02 + u*( 5.0448392644338313e-02 + u*( -2.9713546746201391e-01 + u*( 3.5090739050029035e-01 + u*( -2.5865160758743239e-01 + u*( 1.5034181100350819e-01 + u*( -6.8791765024151963e-02 + u*( 2.3567764847366883e-02 + u*( -5.8664042693176779e-03 + u*( 1.0417420187114196e-03 + u*( -1.2874435074955909e-04 + u*( 1.0549286265432096e-05 + u*( -5.1669973544973548e-07 + u*( 1.1482216343327455e-08))))))))))))));
                     case 12; w = -2.8209730822647999e-01 + u*( 2.5209034433344981e+00 + u*( -4.7643956428425067e+00 + u*( 5.2055434302372445e+00 + u*( -3.9726260291205535e+00 + u*( 2.2119389179101900e+00 + u*( -9.0848269992404529e-01 + u*( 2.7694685078892672e-01 + u*( -6.2866889105902782e-02 + u*( 1.0597339340828923e-02 + u*( -1.3102213541666663e-03 + u*( 1.1561548019881353e-04 + u*( -6.9049873737373750e-06 + u*( 2.5052108385441723e-07 + u*( -4.1753513975736200e-09))))))))))))));
                     case 13; w = 7.0053971489285951e+00 + u*( -2.0151301534481298e+01 + u*( 2.7984344880669187e+01 + u*( -2.3904448146217618e+01 + u*( 1.3816813267601864e+01 + u*( -5.6944785472997719e+00 + u*( 1.7269897884792755e+00 + u*( -3.9237949547223411e-01 + u*( 6.7279900444878471e-02 + u*( -8.6836665185460755e-03 + u*( 8.3211263020833341e-04 + u*( -5.7502417528459207e-05 + u*( 2.7126736111111116e-06 + u*( -7.8287838704505374e-08 + u*( 1.0438378493934050e-09))))))))))))));
                     case 14; w = -2.0910950619651942e+01 + u*( 5.0908492785541867e+01 + u*( -5.5995412042994573e+01 + u*( 3.7171738707356035e+01 + u*( -1.6721280159184964e+01 + u*( 5.4102826988045321e+00 + u*( -1.3015814604582610e+00 + u*( 2.3693401080050075e-01 + u*( -3.2838157371238429e-02 + u*( 3.4518556409832467e-03 + u*( -2.7111665702160499e-04 + u*( 1.5438361792528462e-05 + u*( -6.0281635802469145e-07 + u*( 1.4453139453139453e-08 + u*( -1.6059043836821613e-10))))))))))))));
                     case 15; w = 2.0438514873693794e+01 + u*( -3.8151894430895062e+01 + u*( 3.3064975173442399e+01 + u*( -1.7634653425835946e+01 + u*( 6.4660395894731808e+00 + u*( -1.7242772238595152e+00 + u*( 3.4485544477190305e-01 + u*( -5.2549401108099490e-02 + u*( 6.1307634626116070e-03 + u*( -5.4495675223214287e-04 + u*( 3.6330450148809523e-05 + u*( -1.7614763708513713e-06 + u*( 5.8715879028379031e-08 + u*( -1.2044282877616211e-09 + u*( 1.1470745597729726e-11))))))))))))));
                     otherwise; w = 0;
                end 
                case 16
                switch j
                     case 1; w = 2.6906065416456951e+01 + u*( 5.0448872655856789e+01 + u*( 4.4142763573874696e+01 + u*( 2.3910663602515456e+01 + u*( 8.9664988509432977e+00 + u*( 2.4657871840094066e+00 + u*( 5.1370566333529311e-01 + u*( 8.2559838750314959e-02 + u*( 1.0319979843789366e-02 + u*( 1.0033313737017441e-03 + u*( 7.5249853027630802e-05 + u*( 4.2755598311153861e-06 + u*( 1.7814832629647451e-07 + u*( 5.1388940277829171e-09 + u*( 9.1765964781837809e-11 + u*( 7.6471637318198174e-13)))))))))))))));
                     case 2; w = -3.1182542893014780e+01 + u*( -7.4026716578725555e+01 + u*( -8.0332825660707613e+01 + u*( -5.3145653542702163e+01 + u*( -2.4057637068435678e+01 + u*( -7.9132269620811302e+00 + u*( -1.9574881809719775e+00 + u*( -3.7133290816326520e-01 + u*( -5.4521841143864963e-02 + u*( -6.2013154027042933e-03 + u*( -5.4229129923574369e-04 + u*( -3.5824514991181672e-05 + u*( -1.7313790461938611e-06 + u*( -5.7812557812557817e-08 + u*( -1.1929575421638912e-09 + u*( -1.1470745597729725e-11)))))))))))))));
                     case 3; w = 1.1964435842535368e+01 + u*( 3.3840730260149897e+01 + u*( 4.5512528984647020e+01 + u*( 3.7742658145609525e+01 + u*( 2.1386518775720159e+01 + u*( 8.7496301807760126e+00 + u*( 2.6710832475994524e+00 + u*( 6.2050382653061231e-01 + u*( 1.1078428130511463e-01 + u*( 1.5227256025867140e-02 + u*( 1.6005658436213993e-03 + u*( 1.2651314734648072e-04 + u*( 7.2873799725651587e-06 + u*( 2.8906278906278910e-07 + u*( 7.0659792882015104e-09 + u*( 8.0295219184108080e-11)))))))))))))));
                     case 4; w = -1.1044474882113775e+00 + u*( -5.3659197320903802e+00 + u*( -9.3767810044893576e+00 + u*( -9.8280771783086802e+00 + u*( -7.1559224186307544e+00 + u*( -3.8090439447383933e+00 + u*( -1.5151414609053504e+00 + u*( -4.5595395565633684e-01 + u*( -1.0450727513227519e-01 + u*( -1.8262541642171275e-02 + u*( -2.4182098765432101e-03 + u*( -2.3883009994121106e-04 + u*( -1.7068836513280960e-05 + u*( -8.3507027951472417e-07 + u*( -2.5052108385441724e-08 + u*( -3.4794594979780167e-10)))))))))))))));
                     case 5; w = 3.8996898694583809e-01 + u*( 2.3814204974917877e-01 + u*( 4.3032711372989035e-01 + u*( 7.9628994976215373e-01 + u*( 8.1235292742236864e-01 + u*( 5.7350749559082537e-01 + u*( 3.1092163923182398e-01 + u*( 1.3099489795918348e-01 + u*( 4.2229938271604921e-02 + u*( 1.0269694297472071e-02 + u*( 1.8616255144032923e-03 + u*( 2.4751483084816414e-04 + u*( 2.3459907719166976e-05 + u*( 1.5031265031265033e-06 + u*( 5.8454919566030681e-08 + u*( 1.0438378493934048e-09)))))))))))))));
                     case 6; w = 3.4203960382895132e-01 + u*( -1.5048658352430416e-03 + u*( -1.2884902263375478e-01 + u*( -1.1408913874203518e-02 + u*( 4.6540637860032425e-03 + u*( -1.8805004409171940e-02 + u*( -1.8140860768176267e-02 + u*( -1.0031887755101991e-02 + u*( -4.7789902998236454e-03 + u*( -1.9178057025279254e-03 + u*( -5.7587448559670758e-04 + u*( -1.2180335097001766e-04 + u*( -1.7575445816186560e-05 + u*( -1.6534391534391531e-06 + u*( -9.1857730746619614e-08 + u*( -2.2964432686654912e-09)))))))))))))));
                     case 7; w = 3.4224027010369296e-01 + u*( 1.3122531894081908e-07 + u*( -1.2358153292181553e-01 + u*( 3.9805016588202067e-06 + u*( 2.1773405349795251e-02 + u*( 2.6271310992766949e-05 + u*( -2.4481310013714274e-03 + u*( 5.6295666414727000e-05 + u*( 2.6510141093475493e-04 + u*( 4.3785518322558256e-05 + u*( 1.2602880658436263e-05 + u*( 1.1941504997060561e-05 + u*( 4.7153635116598099e-06 + u*( 9.1857730746619649e-07 + u*( 9.1857730746619628e-08 + u*( 3.8274054477758184e-09)))))))))))))));
                     case 8; w = 3.4224026135533708e-01 + u*( -3.7574792397846952e-15 + u*( -1.2358245149912785e-01 + u*( -5.2197486784961478e-15 + u*( 2.1761463844795555e-02 + u*( -1.3132676641037372e-15 + u*( -2.4919165196942120e-03 + u*( -8.5961082183239520e-17 + u*( 2.0880574452003062e-04 + u*( -6.5503881254332563e-19 + u*( -1.3668430335097271e-05 + u*( 4.2916335994217882e-21 + u*( 7.3486184597295649e-07 + u*( -1.1293772630057338e-22 + u*( -3.9367598891408398e-08 + u*( -4.9209498614260522e-09)))))))))))))));
                     case 9; w = 3.4224026135533925e-01 + u*( -3.4790566415114922e-15 + u*( -1.2358245149911901e-01 + u*( -3.4304886438678900e-15 + u*( 2.1761463844799270e-02 + u*( -1.5441670385216387e-15 + u*( -2.4919165196938512e-03 + u*( -1.5684418734126838e-16 + u*( 2.0880574452004425e-04 + u*( -3.7468220077478220e-18 + u*( -1.3668430335097032e-05 + u*( -1.5811281682080273e-21 + u*( 7.3486184597295564e-07 + u*( 5.6468863150286688e-23 + u*( -3.9367598891408398e-08 + u*( 4.9209498614260522e-09)))))))))))))));
                     case 10; w = 3.4224027010369068e-01 + u*( -1.3122531743168952e-07 + u*( -1.2358153292181757e-01 + u*( -3.9805016563647243e-06 + u*( 2.1773405349794918e-02 + u*( -2.6271310992818977e-05 + u*( -2.4481310013715358e-03 + u*( -5.6295666414787017e-05 + u*( 2.6510141093475314e-04 + u*( -4.3785518322559218e-05 + u*( 1.2602880658436318e-05 + u*( -1.1941504997060561e-05 + u*( 4.7153635116598074e-06 + u*( -9.1857730746619628e-07 + u*( 9.1857730746619614e-08 + u*( -3.8274054477758184e-09)))))))))))))));
                     case 11; w = 3.4203960382895077e-01 + u*( 1.5048658352326478e-03 + u*( -1.2884902263375519e-01 + u*( 1.1408913874202814e-02 + u*( 4.6540637860059140e-03 + u*( 1.8805004409173296e-02 + u*( -1.8140860768175740e-02 + u*( 1.0031887755102066e-02 + u*( -4.7789902998236367e-03 + u*( 1.9178057025279254e-03 + u*( -5.7587448559670768e-04 + u*( 1.2180335097001761e-04 + u*( -1.7575445816186560e-05 + u*( 1.6534391534391533e-06 + u*( -9.1857730746619628e-08 + u*( 2.2964432686654912e-09)))))))))))))));
                     case 12; w = 3.8996898694583926e-01 + u*( -2.3814204974918651e-01 + u*( 4.3032711372988719e-01 + u*( -7.9628994976216305e-01 + u*( 8.1235292742236453e-01 + u*( -5.7350749559082748e-01 + u*( 3.1092163923182370e-01 + u*( -1.3099489795918354e-01 + u*( 4.2229938271604907e-02 + u*( -1.0269694297472069e-02 + u*( 1.8616255144032923e-03 + u*( -2.4751483084816420e-04 + u*( 2.3459907719166972e-05 + u*( -1.5031265031265031e-06 + u*( 5.8454919566030681e-08 + u*( -1.0438378493934048e-09)))))))))))))));
                     case 13; w = -1.1044474882113784e+00 + u*( 5.3659197320903642e+00 + u*( -9.3767810044893540e+00 + u*( 9.8280771783086678e+00 + u*( -7.1559224186307508e+00 + u*( 3.8090439447383884e+00 + u*( -1.5151414609053502e+00 + u*( 4.5595395565633667e-01 + u*( -1.0450727513227513e-01 + u*( 1.8262541642171271e-02 + u*( -2.4182098765432097e-03 + u*( 2.3883009994121103e-04 + u*( -1.7068836513280960e-05 + u*( 8.3507027951472406e-07 + u*( -2.5052108385441720e-08 + u*( 3.4794594979780167e-10)))))))))))))));
                     case 14; w = 1.1964435842535371e+01 + u*( -3.3840730260149904e+01 + u*( 4.5512528984647041e+01 + u*( -3.7742658145609525e+01 + u*( 2.1386518775720170e+01 + u*( -8.7496301807760162e+00 + u*( 2.6710832475994528e+00 + u*( -6.2050382653061242e-01 + u*( 1.1078428130511463e-01 + u*( -1.5227256025867138e-02 + u*( 1.6005658436213993e-03 + u*( -1.2651314734648072e-04 + u*( 7.2873799725651587e-06 + u*( -2.8906278906278910e-07 + u*( 7.0659792882015104e-09 + u*( -8.0295219184108080e-11)))))))))))))));
                     case 15; w = -3.1182542893014780e+01 + u*( 7.4026716578725555e+01 + u*( -8.0332825660707613e+01 + u*( 5.3145653542702171e+01 + u*( -2.4057637068435675e+01 + u*( 7.9132269620811302e+00 + u*( -1.9574881809719775e+00 + u*( 3.7133290816326520e-01 + u*( -5.4521841143864956e-02 + u*( 6.2013154027042933e-03 + u*( -5.4229129923574380e-04 + u*( 3.5824514991181672e-05 + u*( -1.7313790461938611e-06 + u*( 5.7812557812557817e-08 + u*( -1.1929575421638912e-09 + u*( 1.1470745597729725e-11)))))))))))))));
                     case 16; w = 2.6906065416456951e+01 + u*( -5.0448872655856789e+01 + u*( 4.4142763573874696e+01 + u*( -2.3910663602515456e+01 + u*( 8.9664988509432977e+00 + u*( -2.4657871840094066e+00 + u*( 5.1370566333529311e-01 + u*( -8.2559838750314959e-02 + u*( 1.0319979843789366e-02 + u*( -1.0033313737017441e-03 + u*( 7.5249853027630802e-05 + u*( -4.2755598311153861e-06 + u*( 1.7814832629647451e-07 + u*( -5.1388940277829171e-09 + u*( 9.1765964781837809e-11 + u*( -7.6471637318198174e-13)))))))))))))));
                     otherwise; w = 0;
                end 
                case 17
                switch j
                     case 1; w = 3.5488138357040867e+01 + u*( 6.6801201613253383e+01 + u*( 5.8942236717576534e+01 + u*( 3.2360443688081240e+01 + u*( 1.2373110821913412e+01 + u*( 3.4935842320696686e+00 + u*( 7.5351816770130131e-01 + u*( 1.2664170885736156e-01 + u*( 1.6761402642886083e-02 + u*( 1.7528264201710937e-03 + u*( 1.4435041107291360e-04 + u*( 9.2631279832885747e-06 + u*( 4.5407490114159673e-07 + u*( 1.6437100493813461e-08 + u*( 4.1438068471798634e-10 + u*( 6.5000891720468448e-12 + u*( 4.7794773323873859e-14))))))))))))))));
                     case 2; w = -4.5946569342832859e+01 + u*( -1.0692617481314386e+02 + u*( -1.1478513970882069e+02 + u*( -7.5736590532788185e+01 + u*( -3.4468937340463320e+01 + u*( -1.1495871179890882e+01 + u*( -2.9105709330001668e+00 + u*( -5.7128002460958482e-01 + u*( -8.7926857377155873e-02 + u*( -1.0654671063685730e-02 + u*( -1.0136826874203903e-03 + u*( -7.4957460998042604e-05 + u*( -4.2248467089323582e-06 + u*( -1.7551865786819493e-07 + u*( -5.0700695541965378e-09 + u*( -9.1001248408655818e-11 + u*( -7.6471637318198164e-13))))))))))))))));
                     case 3; w = 2.0051771969601770e+01 + u*( 5.5531280725156691e+01 + u*( 7.2665770527679925e+01 + u*( 5.8843550149827671e+01 + u*( 3.2821133000844604e+01 + u*( 1.3349693253822807e+01 + u*( 4.0971523688165155e+00 + u*( 9.6887894282265352e-01 + u*( 1.7863911775534683e-01 + u*( 2.5798795621100976e-02 + u*( 2.9120752632489475e-03 + u*( 2.5447677262455912e-04 + u*( 1.6892732369439548e-05 + u*( 8.2413005590088936e-07 + u*( 2.7885382548080964e-08 + u*( 5.8500802548421594e-10 + u*( 5.7353727988648631e-12))))))))))))))));
                     case 4; w = -2.7355007922911638e+00 + u*( -1.0758967309440882e+01 + u*( -1.7730022246771355e+01 + u*( -1.7855910386070359e+01 + u*( -1.2501275497640611e+01 + u*( -6.4273577273343765e+00 + u*( -2.4951979582358792e+00 + u*( -7.4341984342472223e-01 + u*( -1.7160381579525263e-01 + u*( -3.0806122932531269e-02 + u*( -4.2921870981224288e-03 + u*( -4.5999552767673427e-04 + u*( -3.7233957047325105e-05 + u*( -2.2039364848971332e-06 + u*( -9.0091235924569259e-08 + u*( -2.2750312102163958e-09 + u*( -2.6765073061369358e-11))))))))))))))));
                     case 5; w = 4.8080685754940666e-01 + u*( 6.7679322332555569e-01 + u*( 1.3295786411727306e+00 + u*( 1.9096016458716030e+00 + u*( 1.7738165254286014e+00 + u*( 1.1860246849692007e+00 + u*( 6.0655043196187497e-01 + u*( 2.4126218520948542e-01 + u*( 7.4566691363299284e-02 + u*( 1.7820150086441951e-02 + u*( 3.2718998159400720e-03 + u*( 4.5686349220962934e-04 + u*( 4.7660396645856704e-05 + u*( 3.6008056479699998e-06 + u*( 1.8632505611672282e-07 + u*( 5.9150811465626284e-09 + u*( 8.6986487449450405e-11))))))))))))))));
                     case 6; w = 3.3083166753378568e-01 + u*( -8.8076453172875029e-03 + u*( -1.3956607734760787e-01 + u*( -4.9257978822184585e-02 + u*( -4.5124554644205866e-02 + u*( -6.1249198509294903e-02 + u*( -4.6783507003051804e-02 + u*( -2.5404728653750049e-02 + u*( -1.1147673807026425e-02 + u*( -3.9485775758630074e-03 + u*( -1.0818457165209187e-03 + u*( -2.2164230505701866e-04 + u*( -3.3114103028744229e-05 + u*( -3.5002492684784338e-06 + u*( -2.4843340815563032e-07 + u*( -1.0647146063812730e-08 + u*( -2.0876756987868101e-10))))))))))))))));
                     case 7; w = 3.3220887968440743e-01 + u*( 6.5124466913968389e-06 + u*( -1.1312360405572261e-01 + u*( 1.0130465602479965e-04 + u*( 1.9042513877436530e-02 + u*( 3.5118727148732571e-04 + u*( -1.6098907638155309e-03 + u*( 4.0876634009962768e-04 + u*( 4.6839894020579372e-04 + u*( 1.8158162315289693e-04 + u*( 7.4598859203532737e-05 + u*( 3.0672875101043510e-05 + u*( 8.9384269975994522e-06 + u*( 1.6754467347638653e-06 + u*( 1.9519767783656676e-07 + u*( 1.3013178522437781e-08 + u*( 3.8274054477758185e-10))))))))))))))));
                     case 8; w = 3.3220826914247464e-01 + u*( -5.8414547797179246e-13 + u*( -1.1315616629211941e-01 + u*( -7.9440383243409664e-11 + u*( 1.8823020283950668e-02 + u*( -2.4780977890441289e-09 + u*( -2.0391226799702854e-03 + u*( -2.5961000765497433e-08 + u*( 1.6180471438067897e-04 + u*( -1.0384400271793464e-07 + u*( -1.0187692135753885e-05 + u*( -1.5859811324222923e-07 + u*( 3.7412888252008622e-07 + u*( -8.1332365765236104e-08 + u*( -5.5770765096161928e-08 + u*( -9.2951275160269858e-09 + u*( -5.4677220682511691e-10))))))))))))))));
                     case 9; w = 3.3220826914249191e-01 + u*( -2.0821295796458210e-15 + u*( -1.1315616628360614e-01 + u*( -4.3582649956798410e-15 + u*( 1.8823020800222107e-02 + u*( -1.7796214631976461e-15 + u*( -2.0391135936198559e-03 + u*( -1.4649688932647301e-16 + u*( 1.6186312663221314e-04 + u*( -1.9157767684524138e-18 + u*( -1.0042310531948603e-05 + u*( 1.3340768919255230e-20 + u*( 5.0629397688859398e-07 + u*( -1.0587911840678754e-22 + u*( -2.0914036911060677e-08 + u*( 0.0000000000000000e+00 + u*( 6.1511873267825652e-10))))))))))))))));
                     case 10; w = 3.3220826914247864e-01 + u*( 5.7053311700378333e-13 + u*( -1.1315616629210053e-01 + u*( 7.9426025676613272e-11 + u*( 1.8823020283958776e-02 + u*( 2.4780936643524869e-09 + u*( -2.0391226799694570e-03 + u*( 2.5961000438699311e-08 + u*( 1.6180471438071180e-04 + u*( 1.0384400271199948e-07 + u*( -1.0187692135753471e-05 + u*( 1.5859811324221483e-07 + u*( 3.7412888252008405e-07 + u*( 8.1332365765236051e-08 + u*( -5.5770765096161928e-08 + u*( 9.2951275160269842e-09 + u*( -5.4677220682511691e-10))))))))))))))));
                     case 11; w = 3.3220887968440238e-01 + u*( -6.5124466935883334e-06 + u*( -1.1312360405572433e-01 + u*( -1.0130465601738547e-04 + u*( 1.9042513877441204e-02 + u*( -3.5118727148499907e-04 + u*( -1.6098907638148162e-03 + u*( -4.0876634009955650e-04 + u*( 4.6839894020580928e-04 + u*( -1.8158162315289755e-04 + u*( 7.4598859203532601e-05 + u*( -3.0672875101043516e-05 + u*( 8.9384269975994455e-06 + u*( -1.6754467347638645e-06 + u*( 1.9519767783656674e-07 + u*( -1.3013178522437780e-08 + u*( 3.8274054477758185e-10))))))))))))))));
                     case 12; w = 3.3083166753378562e-01 + u*( 8.8076453172697845e-03 + u*( -1.3956607734762064e-01 + u*( 4.9257978822168154e-02 + u*( -4.5124554644212340e-02 + u*( 6.1249198509292536e-02 + u*( -4.6783507003052081e-02 + u*( 2.5404728653749962e-02 + u*( -1.1147673807026427e-02 + u*( 3.9485775758630074e-03 + u*( -1.0818457165209185e-03 + u*( 2.2164230505701858e-04 + u*( -3.3114103028744229e-05 + u*( 3.5002492684784342e-06 + u*( -2.4843340815563038e-07 + u*( 1.0647146063812729e-08 + u*( -2.0876756987868101e-10))))))))))))))));
                     case 13; w = 4.8080685754940455e-01 + u*( -6.7679322332557401e-01 + u*( 1.3295786411727237e+00 + u*( -1.9096016458716181e+00 + u*( 1.7738165254285965e+00 + u*( -1.1860246849692049e+00 + u*( 6.0655043196187453e-01 + u*( -2.4126218520948550e-01 + u*( 7.4566691363299298e-02 + u*( -1.7820150086441951e-02 + u*( 3.2718998159400716e-03 + u*( -4.5686349220962934e-04 + u*( 4.7660396645856704e-05 + u*( -3.6008056479699990e-06 + u*( 1.8632505611672277e-07 + u*( -5.9150811465626284e-09 + u*( 8.6986487449450405e-11))))))))))))))));
                     case 14; w = -2.7355007922911554e+00 + u*( 1.0758967309440871e+01 + u*( -1.7730022246771309e+01 + u*( 1.7855910386070356e+01 + u*( -1.2501275497640599e+01 + u*( 6.4273577273343703e+00 + u*( -2.4951979582358779e+00 + u*( 7.4341984342472178e-01 + u*( -1.7160381579525252e-01 + u*( 3.0806122932531266e-02 + u*( -4.2921870981224271e-03 + u*( 4.5999552767673422e-04 + u*( -3.7233957047325098e-05 + u*( 2.2039364848971332e-06 + u*( -9.0091235924569246e-08 + u*( 2.2750312102163958e-09 + u*( -2.6765073061369358e-11))))))))))))))));
                     case 15; w = 2.0051771969601781e+01 + u*( -5.5531280725156705e+01 + u*( 7.2665770527679967e+01 + u*( -5.8843550149827678e+01 + u*( 3.2821133000844625e+01 + u*( -1.3349693253822810e+01 + u*( 4.0971523688165172e+00 + u*( -9.6887894282265352e-01 + u*( 1.7863911775534685e-01 + u*( -2.5798795621100976e-02 + u*( 2.9120752632489470e-03 + u*( -2.5447677262455912e-04 + u*( 1.6892732369439548e-05 + u*( -8.2413005590088936e-07 + u*( 2.7885382548080967e-08 + u*( -5.8500802548421594e-10 + u*( 5.7353727988648631e-12))))))))))))))));
                     case 16; w = -4.5946569342832859e+01 + u*( 1.0692617481314386e+02 + u*( -1.1478513970882071e+02 + u*( 7.5736590532788171e+01 + u*( -3.4468937340463320e+01 + u*( 1.1495871179890882e+01 + u*( -2.9105709330001672e+00 + u*( 5.7128002460958482e-01 + u*( -8.7926857377155859e-02 + u*( 1.0654671063685731e-02 + u*( -1.0136826874203905e-03 + u*( 7.4957460998042604e-05 + u*( -4.2248467089323573e-06 + u*( 1.7551865786819493e-07 + u*( -5.0700695541965378e-09 + u*( 9.1001248408655818e-11 + u*( -7.6471637318198164e-13))))))))))))))));
                     case 17; w = 3.5488138357040867e+01 + u*( -6.6801201613253383e+01 + u*( 5.8942236717576534e+01 + u*( -3.2360443688081240e+01 + u*( 1.2373110821913412e+01 + u*( -3.4935842320696686e+00 + u*( 7.5351816770130131e-01 + u*( -1.2664170885736156e-01 + u*( 1.6761402642886083e-02 + u*( -1.7528264201710937e-03 + u*( 1.4435041107291360e-04 + u*( -9.2631279832885747e-06 + u*( 4.5407490114159673e-07 + u*( -1.6437100493813461e-08 + u*( 4.1438068471798634e-10 + u*( -6.5000891720468448e-12 + u*( 4.7794773323873859e-14))))))))))))))));
                     otherwise; w = 0;
                end 
                case 18
                switch j
                     case 1; w = 4.6887183471565933e+01 + u*( 8.8564679890735661e+01 + u*( 7.8724159902876181e+01 + u*( 4.3735644390486755e+01 + u*( 1.7008306151855958e+01 + u*( 4.9135106660917209e+00 + u*( 1.0918912591314938e+00 + u*( 1.9064768016581635e-01 + u*( 2.6478844467474488e-02 + u*( 2.9420938297193872e-03 + u*( 2.6151945153061220e-04 + u*( 1.8491274350649356e-05 + u*( 1.0272930194805198e-06 + u*( 4.3901411088911095e-08 + u*( 1.3936955901241619e-09 + u*( 3.0971013113870259e-11 + u*( 4.3015295991486473e-13 + u*( 2.8114572543455210e-15)))))))))))))))));
                     case 2; w = -6.7067917115781142e+01 + u*( -1.5358990885737700e+02 + u*( -1.6343042884523643e+02 + u*( -1.0761097357708366e+02 + u*( -4.9205839208956064e+01 + u*( -1.6606086576172189e+01 + u*( -4.2880080514344820e+00 + u*( -8.6611825583821489e-01 + u*( -1.3864083303315541e-01 + u*( -1.7697865857859346e-02 + u*( -1.8024765172272616e-03 + u*( -1.4569022316418152e-04 + u*( -9.2340505751964102e-06 + u*( -4.4943241557824901e-07 + u*( -1.6225369647988698e-08 + u*( -4.0950561783895120e-10 + u*( -6.4522943987229703e-12 + u*( -4.7794773323873853e-14)))))))))))))))));
                     case 3; w = 3.2998787042332253e+01 + u*( 8.9429229812327023e+01 + u*( 1.1430572963442528e+02 + u*( 9.0771996765531881e+01 + u*( 4.9985645962351690e+01 + u*( 2.0236465058884974e+01 + u*( 6.2384352728675649e+00 + u*( 1.4969608577806119e+00 + u*( 2.8333758011306381e-01 + u*( 4.2584764591600516e-02 + u*( 5.0869669627110117e-03 + u*( 4.8062282046657049e-04 + u*( 3.5502595398428734e-05 + u*( 2.0086250555000557e-06 + u*( 8.4103506722554340e-08 + u*( 2.4570337070337067e-09 + u*( 4.4735907831145924e-11 + u*( 3.8235818659099087e-13)))))))))))))))));
                     case 4; w = -5.8334938196628920e+00 + u*( -2.0595565963325868e+01 + u*( -3.2393998066445292e+01 + u*( -3.1477776318526942e+01 + u*( -2.1326721670015953e+01 + u*( -1.0665560915141000e+01 + u*( -4.0622400518077635e+00 + u*( -1.2008350605867346e+00 + u*( -2.7870323621346699e-01 + u*( -5.1088704796154585e-02 + u*( -7.4028289556563372e-03 + u*( -8.4405250420875426e-04 + u*( -7.4887014991181666e-05 + u*( -5.0676320207570222e-06 + u*( -2.5286111595635406e-07 + u*( -8.7751203822632388e-09 + u*( -1.8926730236254048e-10 + u*( -1.9117909329549544e-12)))))))))))))))));
                     case 5; w = 7.3011946028804298e-01 + u*( 1.7207191885073605e+00 + u*( 3.3120581764877985e+00 + u*( 4.2282799244061851e+00 + u*( 3.6675177000372154e+00 + u*( 2.3314435572866459e+00 + u*( 1.1365617371632983e+00 + u*( 4.3307407308988366e-01 + u*( 1.2977404720568775e-01 + u*( 3.0606751887676358e-02 + u*( 5.6684441137566143e-03 + u*( 8.1956406826198499e-04 + u*( 9.1474642255892253e-05 + u*( 7.7294185367102045e-06 + u*( 4.7839891589891597e-07 + u*( 2.0475280891947560e-08 + u*( 5.4199272949272949e-10 + u*( 6.6912682653423394e-12)))))))))))))))));
                     case 6; w = 3.1628105178297428e-01 + u*( -3.8094047639235115e-02 + u*( -2.0556829580528790e-01 + u*( -1.6875316596024745e-01 + u*( -1.7988625403339287e-01 + u*( -1.6936901285925418e-01 + u*( -1.1384454790965244e-01 + u*( -5.8156967474490014e-02 + u*( -2.3735652970679037e-02 + u*( -7.7706731564153439e-03 + u*( -2.0070408950617282e-03 + u*( -4.0153581950456947e-04 + u*( -6.1162843714927047e-05 + u*( -6.9472628066378051e-06 + u*( -5.6993546576879891e-07 + u*( -3.1941438191438190e-08 + u*( -1.0960297418630753e-09 + u*( -1.7397297489890085e-11)))))))))))))))));
                     case 7; w = 3.2302112128379673e-01 + u*( 9.9679532040898404e-05 + u*( -1.0371835668188797e-01 + u*( 9.9673257873929798e-04 + u*( 1.8155294262070173e-02 + u*( 2.2669956634782254e-03 + u*( 5.7945777216681944e-04 + u*( 1.7794164540816664e-03 + u*( 1.2378403328924701e-03 + u*( 5.5382461144179852e-04 + u*( 2.1282517636684451e-04 + u*( 6.9344862313612285e-05 + u*( 1.7317269921436599e-05 + u*( 3.1142902236652256e-06 + u*( 3.8830767997434672e-07 + u*( 3.1941438191438196e-08 + u*( 1.5657567740901075e-09 + u*( 3.4794594979780171e-11)))))))))))))))));
                     case 8; w = 3.2300939403396972e-01 + u*( -2.0914320234219767e-09 + u*( -1.0411708317588460e-01 + u*( -8.3656157400054738e-08 + u*( 1.6410865850960656e-02 + u*( -7.6127094311473153e-07 + u*( -1.6882991622584832e-03 + u*( -2.3925658225746768e-06 + u*( 1.2420969545223195e-04 + u*( -2.9907072782835146e-06 + u*( -9.9009511211893845e-06 + u*( -1.5225418871252779e-06 + u*( -3.9958112874779922e-07 + u*( -2.9279651675484939e-07 + u*( -9.8418997228521077e-08 + u*( -1.6731229528848580e-08 + u*( -1.4762849584278154e-09 + u*( -5.4677220682511690e-11)))))))))))))))));
                     case 9; w = 3.2300939415699759e-01 + u*( -1.6557887785903664e-14 + u*( -1.0411706644463017e-01 + u*( -1.5194906225872946e-14 + u*( 1.6411158647488929e-02 + u*( -2.9949034354374940e-15 + u*( -1.6867766203703292e-03 + u*( -1.5649026189047218e-16 + u*( 1.2720040273051054e-04 + u*( -1.3903983547344982e-18 + u*( -7.5083852985646201e-06 + u*( -1.2400313220465530e-20 + u*( 3.6168981481481591e-07 + u*( 1.8684550307080154e-22 + u*( -1.4762849584278137e-08 + u*( 3.8926146473083657e-24 + u*( 6.1511873267825632e-10 + u*( 6.8346525853139614e-11)))))))))))))))));
                     case 10; w = 3.2300939415699981e-01 + u*( 8.7389246948726411e-15 + u*( -1.0411706644461764e-01 + u*( 1.0448431134664892e-15 + u*( 1.6411158647494768e-02 + u*( -2.0212130567887399e-15 + u*( -1.6867766203697257e-03 + u*( -2.1635174319791092e-16 + u*( 1.2720040273052978e-04 + u*( -3.0000225382219679e-18 + u*( -7.5083852985645735e-06 + u*( 3.5021075458903902e-20 + u*( 3.6168981481481199e-07 + u*( -5.6053650921240462e-22 + u*( -1.4762849584278137e-08 + u*( -5.4496605062317115e-24 + u*( 6.1511873267825632e-10 + u*( -6.8346525853139614e-11)))))))))))))))));
                     case 11; w = 3.2300939403397638e-01 + u*( 2.0914073809532974e-09 + u*( -1.0411708317583972e-01 + u*( 8.3656144963930003e-08 + u*( 1.6410865850984082e-02 + u*( 7.6127094061522582e-07 + u*( -1.6882991622559212e-03 + u*( 2.3925658222958573e-06 + u*( 1.2420969545230977e-04 + u*( 2.9907072782752954e-06 + u*( -9.9009511211888949e-06 + u*( 1.5225418871252180e-06 + u*( -3.9958112874780224e-07 + u*( 2.9279651675484981e-07 + u*( -9.8418997228521011e-08 + u*( 1.6731229528848577e-08 + u*( -1.4762849584278154e-09 + u*( 5.4677220682511690e-11)))))))))))))))));
                     case 12; w = 3.2302112128378618e-01 + u*( -9.9679532064481156e-05 + u*( -1.0371835668192195e-01 + u*( -9.9673257875976013e-04 + u*( 1.8155294262061295e-02 + u*( -2.2669956634808245e-03 + u*( 5.7945777216657245e-04 + u*( -1.7794164540817813e-03 + u*( 1.2378403328924734e-03 + u*( -5.5382461144179917e-04 + u*( 2.1282517636684432e-04 + u*( -6.9344862313612312e-05 + u*( 1.7317269921436585e-05 + u*( -3.1142902236652243e-06 + u*( 3.8830767997434651e-07 + u*( -3.1941438191438196e-08 + u*( 1.5657567740901071e-09 + u*( -3.4794594979780171e-11)))))))))))))))));
                     case 13; w = 3.1628105178297006e-01 + u*( 3.8094047639204556e-02 + u*( -2.0556829580530958e-01 + u*( 1.6875316596022064e-01 + u*( -1.7988625403340411e-01 + u*( 1.6936901285924907e-01 + u*( -1.1384454790965300e-01 + u*( 5.8156967474489916e-02 + u*( -2.3735652970679044e-02 + u*( 7.7706731564153421e-03 + u*( -2.0070408950617282e-03 + u*( 4.0153581950456936e-04 + u*( -6.1162843714927047e-05 + u*( 6.9472628066378068e-06 + u*( -5.6993546576879891e-07 + u*( 3.1941438191438190e-08 + u*( -1.0960297418630753e-09 + u*( 1.7397297489890085e-11)))))))))))))))));
                     case 14; w = 7.3011946028805430e-01 + u*( -1.7207191885073507e+00 + u*( 3.3120581764878350e+00 + u*( -4.2282799244061859e+00 + u*( 3.6675177000372186e+00 + u*( -2.3314435572866503e+00 + u*( 1.1365617371632981e+00 + u*( -4.3307407308988366e-01 + u*( 1.2977404720568786e-01 + u*( -3.0606751887676344e-02 + u*( 5.6684441137566152e-03 + u*( -8.1956406826198477e-04 + u*( 9.1474642255892240e-05 + u*( -7.7294185367102011e-06 + u*( 4.7839891589891576e-07 + u*( -2.0475280891947560e-08 + u*( 5.4199272949272949e-10 + u*( -6.6912682653423394e-12)))))))))))))))));
                     case 15; w = -5.8334938196628690e+00 + u*( 2.0595565963325857e+01 + u*( -3.2393998066445199e+01 + u*( 3.1477776318526928e+01 + u*( -2.1326721670015928e+01 + u*( 1.0665560915140997e+01 + u*( -4.0622400518077590e+00 + u*( 1.2008350605867342e+00 + u*( -2.7870323621346682e-01 + u*( 5.1088704796154571e-02 + u*( -7.4028289556563372e-03 + u*( 8.4405250420875415e-04 + u*( -7.4887014991181639e-05 + u*( 5.0676320207570205e-06 + u*( -2.5286111595635406e-07 + u*( 8.7751203822632388e-09 + u*( -1.8926730236254048e-10 + u*( 1.9117909329549544e-12)))))))))))))))));
                     case 16; w = 3.2998787042332268e+01 + u*( -8.9429229812327080e+01 + u*( 1.1430572963442532e+02 + u*( -9.0771996765531924e+01 + u*( 4.9985645962351711e+01 + u*( -2.0236465058884985e+01 + u*( 6.2384352728675667e+00 + u*( -1.4969608577806122e+00 + u*( 2.8333758011306381e-01 + u*( -4.2584764591600516e-02 + u*( 5.0869669627110117e-03 + u*( -4.8062282046657038e-04 + u*( 3.5502595398428734e-05 + u*( -2.0086250555000557e-06 + u*( 8.4103506722554340e-08 + u*( -2.4570337070337071e-09 + u*( 4.4735907831145924e-11 + u*( -3.8235818659099087e-13)))))))))))))))));
                     case 17; w = -6.7067917115781142e+01 + u*( 1.5358990885737700e+02 + u*( -1.6343042884523641e+02 + u*( 1.0761097357708366e+02 + u*( -4.9205839208956050e+01 + u*( 1.6606086576172192e+01 + u*( -4.2880080514344820e+00 + u*( 8.6611825583821489e-01 + u*( -1.3864083303315539e-01 + u*( 1.7697865857859346e-02 + u*( -1.8024765172272614e-03 + u*( 1.4569022316418155e-04 + u*( -9.2340505751964102e-06 + u*( 4.4943241557824901e-07 + u*( -1.6225369647988698e-08 + u*( 4.0950561783895115e-10 + u*( -6.4522943987229703e-12 + u*( 4.7794773323873853e-14)))))))))))))))));
                     case 18; w = 4.6887183471565933e+01 + u*( -8.8564679890735661e+01 + u*( 7.8724159902876181e+01 + u*( -4.3735644390486755e+01 + u*( 1.7008306151855958e+01 + u*( -4.9135106660917209e+00 + u*( 1.0918912591314938e+00 + u*( -1.9064768016581635e-01 + u*( 2.6478844467474488e-02 + u*( -2.9420938297193872e-03 + u*( 2.6151945153061220e-04 + u*( -1.8491274350649356e-05 + u*( 1.0272930194805198e-06 + u*( -4.3901411088911095e-08 + u*( 1.3936955901241619e-09 + u*( -3.0971013113870259e-11 + u*( 4.3015295991486473e-13 + u*( -2.8114572543455210e-15)))))))))))))))));
                     otherwise; w = 0;
                end 
                case 19
                switch j
                     case 1; w = 6.2041726508848370e+01 + u*( 1.1755274496413374e+02 + u*( 1.0517877181001440e+02 + u*( 5.9047731542464234e+01 + u*( 2.3308315082551665e+01 + u*( 6.8698191822257550e+00 + u*( 1.5668008661216632e+00 + u*( 2.8273098336030006e-01 + u*( 4.0921589696885537e-02 + u*( 4.7861508417409975e-03 + u*( 4.5342481658598943e-04 + u*( 3.4711947681224070e-05 + u*( 2.1314353839348120e-06 + u*( 1.0355151662841189e-07 + u*( 3.8929141589628539e-09 + u*( 1.0927478340948359e-10 + u*( 2.1567391462398081e-12 + u*( 2.6708843916282449e-14 + u*( 1.5619206968586228e-16))))))))))))))))));
                     case 2; w = -9.7162005287321065e+01 + u*( -2.1958456942775456e+02 + u*( -2.3195854258187370e+02 + u*( -1.5248940689950487e+02 + u*( -7.0016893053611184e+01 + u*( -2.3872602321451406e+01 + u*( -6.2695026544234951e+00 + u*( -1.2976999787664547e+00 + u*( -2.1473636005891306e-01 + u*( -2.8633188995618286e-02 + u*( -3.0850935191344062e-03 + u*( -2.6804897943923762e-04 + u*( -1.8646275300802758e-05 + u*( -1.0246499685157091e-06 + u*( -4.3510509586588465e-08 + u*( -1.3778914517450677e-09 + u*( -3.0648398393934109e-11 + u*( -4.2734150266051919e-13 + u*( -2.8114572543455210e-15))))))))))))))))));
                     case 3; w = 5.3414714183682392e+01 + u*( 1.4179955730265365e+02 + u*( 1.7761013437925556e+02 + u*( 1.3875943005063152e+02 + u*( 7.5607525421457083e+01 + u*( 3.0493847242574052e+01 + u*( 9.4363605529616397e+00 + u*( 2.2922116114930042e+00 + u*( 4.4341409815532135e-01 + u*( 6.8870582591675680e-02 + u*( 8.6153590713408702e-03 + u*( 8.6654036266745602e-04 + u*( 6.9599562418606745e-05 + u*( 4.4058631219094897e-06 + u*( 2.1508535186223060e-07 + u*( 7.8166280664351651e-09 + u*( 1.9921458956057168e-10 + u*( 3.1783524260376107e-12 + u*( 2.3897386661936929e-14))))))))))))))))));
                     case 4; w = -1.1511962885627518e+01 + u*( -3.7997394581589084e+01 + u*( -5.7508956546292701e+01 + u*( -5.4158798401100370e+01 + u*( -3.5691452531465174e+01 + u*( -1.7450327875607844e+01 + u*( -6.5450311530989920e+00 + u*( -1.9226609263691397e+00 + u*( -4.4819355408474759e-01 + u*( -8.3540981893806252e-02 + u*( -1.2487780626648928e-02 + u*( -1.4946500930516828e-03 + u*( -1.4230214514593137e-04 + u*( -1.0640411971430490e-05 + u*( -6.1163305985974634e-07 + u*( -2.6100024722158764e-08 + u*( -7.7915039472579164e-10 + u*( -1.4529611090457651e-11 + u*( -1.2745272886366362e-13))))))))))))))))));
                     case 5; w = 1.3282185064293035e+00 + u*( 4.0250172469605587e+00 + u*( 7.4347708251021203e+00 + u*( 8.8169372317673549e+00 + u*( 7.2465490363991547e+00 + u*( 4.4090183771231057e+00 + u*( 2.0662264616131947e+00 + u*( 7.6136742107362110e-01 + u*( 2.2281353277594215e-01 + u*( 5.2016005350777640e-02 + u*( 9.6942718315557019e-03 + u*( 1.4385138683638060e-03 + u*( 1.6879100227692347e-04 + u*( 1.5465306693424460e-05 + u*( 1.0835434768191465e-06 + u*( 5.6090352813787567e-08 + u*( 2.0227942939996508e-09 + u*( 4.5405034657680168e-11 + u*( 4.7794773323873859e-13))))))))))))))))));
                     case 6; w = 2.8832407930647252e-01 + u*( -1.3456046153066603e-01 + u*( -4.2220929093695125e-01 + u*( -4.9503920205670821e-01 + u*( -5.1343132512090728e-01 + u*( -4.1941384782272317e-01 + u*( -2.5857423928664436e-01 + u*( -1.2427094117393832e-01 + u*( -4.7798189021922947e-02 + u*( -1.4801703735115008e-02 + u*( -3.6692699856228283e-03 + u*( -7.2125046572565367e-04 + u*( -1.1117844843837680e-04 + u*( -1.3249508764555056e-05 + u*( -1.1954101309570057e-06 + u*( -7.8958749869243675e-08 + u*( -3.6042516511266514e-09 + u*( -1.0170727763320357e-10 + u*( -1.3382536530684680e-12))))))))))))))))));
                     case 7; w = 3.1464890998232264e-01 + u*( 8.2438194495385887e-04 + u*( -9.3417528210308395e-02 + u*( 5.9768173361442927e-03 + u*( 2.3371552800010860e-02 + u*( 1.0028454514000400e-02 + u*( 7.2709954932301257e-03 + u*( 5.9389697386525962e-03 + u*( 3.3557045508805085e-03 + u*( 1.4376275578384565e-03 + u*( 5.0655806113663162e-04 + u*( 1.4645406347111732e-04 + u*( 3.3438973094418331e-05 + u*( 5.8209204485607867e-06 + u*( 7.5055203364665197e-07 + u*( 6.9305034100558800e-08 + u*( 4.3384510615413380e-09 + u*( 1.6527432615395580e-10 + u*( 2.8995495816483473e-12))))))))))))))))));
                     case 8; w = 3.1453438341167761e-01 + u*( -2.0936371923447998e-07 + u*( -9.6221138659803807e-02 + u*( -4.2182893990920751e-06 + u*( 1.4399999361681980e-02 + u*( -1.9685336917280484e-05 + u*( -1.4373923775653223e-03 + u*( -3.2496229892493892e-05 + u*( 7.1398268180908548e-05 + u*( -2.2064123361378571e-05 + u*( -1.8930944095302139e-05 + u*( -6.4154653236273053e-06 + u*( -2.2305836243553305e-06 + u*( -7.6422848413589036e-07 + u*( -1.9018352816715919e-07 + u*( -3.1040092492914368e-08 + u*( -3.1874334329691481e-09 + u*( -1.8888494417594947e-10 + u*( -4.9706564256828813e-12))))))))))))))))));
                     case 9; w = 3.1453440085864343e-01 + u*( -1.6539911881522266e-14 + u*( -9.6219952265849465e-02 + u*( 3.5630364741257783e-13 + u*( 1.4410545085956057e-02 + u*( 1.5091403540050443e-11 + u*( -1.4089579802286443e-03 + u*( 2.2420859224632530e-10 + u*( 1.0118668443996983e-04 + u*( 1.3701639196883930e-09 + u*( -5.6916479801290546e-06 + u*( 3.5873382756581210e-09 + u*( 2.6571463305132801e-07 + u*( 3.8632873738500492e-09 + u*( -7.3045349505542910e-09 + u*( 1.4717285233709395e-09 + u*( 8.7654419406651595e-10 + u*( 1.2985839912096526e-10 + u*( 6.8346525853139620e-12))))))))))))))))));
                     case 10; w = 3.1453440085865725e-01 + u*( -2.3689298206042762e-15 + u*( -9.6219952265824443e-02 + u*( -7.8662194542518841e-15 + u*( 1.4410545083283821e-02 + u*( -3.0051645555908032e-15 + u*( -1.4089580456210877e-03 + u*( -2.1039575931868026e-16 + u*( 1.0118606786621452e-04 + u*( -2.0765247988868966e-18 + u*( -5.6941142751941387e-06 + u*( 7.7056469507162044e-22 + u*( 2.6152940506299188e-07 + u*( -5.8233515123733148e-22 + u*( -1.0064025931874875e-08 + u*( -1.4705433112053825e-24 + u*( 3.2464599780241234e-10 + u*( 0.0000000000000000e+00 + u*( -7.5940584281266239e-12))))))))))))))))));
                     case 11; w = 3.1453440085865236e-01 + u*( 9.1808310988993074e-15 + u*( -9.6219952265788125e-02 + u*( -3.5920636742466043e-13 + u*( 1.4410545085985681e-02 + u*( -1.5093714956399658e-11 + u*( -1.4089579802257287e-03 + u*( -2.2420887152558764e-10 + u*( 1.0118668444003938e-04 + u*( -1.3701639250946691e-09 + u*( -5.6916479801291927e-06 + u*( -3.5873382756595446e-09 + u*( 2.6571463305132182e-07 + u*( -3.8632873738497960e-09 + u*( -7.3045349505542025e-09 + u*( -1.4717285233709412e-09 + u*( 8.7654419406651595e-10 + u*( -1.2985839912096526e-10 + u*( 6.8346525853139620e-12))))))))))))))))));
                     case 12; w = 3.1453438341167728e-01 + u*( 2.0936362572978946e-07 + u*( -9.6221138659790983e-02 + u*( 4.2182893263343327e-06 + u*( 1.4399999361696204e-02 + u*( 1.9685336905532310e-05 + u*( -1.4373923775629904e-03 + u*( 3.2496229891812762e-05 + u*( 7.1398268181002399e-05 + u*( 2.2064123361364361e-05 + u*( -1.8930944095301482e-05 + u*( 6.4154653236271571e-06 + u*( -2.2305836243553381e-06 + u*( 7.6422848413588993e-07 + u*( -1.9018352816715929e-07 + u*( 3.1040092492914348e-08 + u*( -3.1874334329691481e-09 + u*( 1.8888494417594944e-10 + u*( -4.9706564256828813e-12))))))))))))))))));
                     case 13; w = 3.1464890998230444e-01 + u*( -8.2438194499792584e-04 + u*( -9.3417528210365100e-02 + u*( -5.9768173361865948e-03 + u*( 2.3371552799992312e-02 + u*( -1.0028454514006923e-02 + u*( 7.2709954932294396e-03 + u*( -5.9389697386527975e-03 + u*( 3.3557045508805015e-03 + u*( -1.4376275578384602e-03 + u*( 5.0655806113663107e-04 + u*( -1.4645406347111727e-04 + u*( 3.3438973094418318e-05 + u*( -5.8209204485607825e-06 + u*( 7.5055203364665176e-07 + u*( -6.9305034100558761e-08 + u*( 4.3384510615413380e-09 + u*( -1.6527432615395578e-10 + u*( 2.8995495816483473e-12))))))))))))))))));
                     case 14; w = 2.8832407930649046e-01 + u*( 1.3456046153067810e-01 + u*( -4.2220929093692555e-01 + u*( 4.9503920205669893e-01 + u*( -5.1343132512091105e-01 + u*( 4.1941384782271868e-01 + u*( -2.5857423928664391e-01 + u*( 1.2427094117393854e-01 + u*( -4.7798189021922816e-02 + u*( 1.4801703735115027e-02 + u*( -3.6692699856228270e-03 + u*( 7.2125046572565378e-04 + u*( -1.1117844843837678e-04 + u*( 1.3249508764555056e-05 + u*( -1.1954101309570057e-06 + u*( 7.8958749869243662e-08 + u*( -3.6042516511266514e-09 + u*( 1.0170727763320357e-10 + u*( -1.3382536530684680e-12))))))))))))))))));
                     case 15; w = 1.3282185064293344e+00 + u*( -4.0250172469605285e+00 + u*( 7.4347708251022038e+00 + u*( -8.8169372317673389e+00 + u*( 7.2465490363991796e+00 + u*( -4.4090183771231022e+00 + u*( 2.0662264616131969e+00 + u*( -7.6136742107362088e-01 + u*( 2.2281353277594226e-01 + u*( -5.2016005350777640e-02 + u*( 9.6942718315556984e-03 + u*( -1.4385138683638060e-03 + u*( 1.6879100227692342e-04 + u*( -1.5465306693424456e-05 + u*( 1.0835434768191463e-06 + u*( -5.6090352813787540e-08 + u*( 2.0227942939996512e-09 + u*( -4.5405034657680168e-11 + u*( 4.7794773323873859e-13))))))))))))))))));
                     case 16; w = -1.1511962885627492e+01 + u*( 3.7997394581588985e+01 + u*( -5.7508956546292630e+01 + u*( 5.4158798401100299e+01 + u*( -3.5691452531465139e+01 + u*( 1.7450327875607833e+01 + u*( -6.5450311530989858e+00 + u*( 1.9226609263691390e+00 + u*( -4.4819355408474748e-01 + u*( 8.3540981893806238e-02 + u*( -1.2487780626648924e-02 + u*( 1.4946500930516824e-03 + u*( -1.4230214514593135e-04 + u*( 1.0640411971430487e-05 + u*( -6.1163305985974634e-07 + u*( 2.6100024722158764e-08 + u*( -7.7915039472579164e-10 + u*( 1.4529611090457651e-11 + u*( -1.2745272886366362e-13))))))))))))))))));
                     case 17; w = 5.3414714183682420e+01 + u*( -1.4179955730265374e+02 + u*( 1.7761013437925570e+02 + u*( -1.3875943005063158e+02 + u*( 7.5607525421457126e+01 + u*( -3.0493847242574059e+01 + u*( 9.4363605529616432e+00 + u*( -2.2922116114930047e+00 + u*( 4.4341409815532129e-01 + u*( -6.8870582591675680e-02 + u*( 8.6153590713408702e-03 + u*( -8.6654036266745592e-04 + u*( 6.9599562418606745e-05 + u*( -4.4058631219094897e-06 + u*( 2.1508535186223060e-07 + u*( -7.8166280664351651e-09 + u*( 1.9921458956057170e-10 + u*( -3.1783524260376111e-12 + u*( 2.3897386661936929e-14))))))))))))))))));
                     case 18; w = -9.7162005287321065e+01 + u*( 2.1958456942775453e+02 + u*( -2.3195854258187370e+02 + u*( 1.5248940689950484e+02 + u*( -7.0016893053611184e+01 + u*( 2.3872602321451403e+01 + u*( -6.2695026544234942e+00 + u*( 1.2976999787664549e+00 + u*( -2.1473636005891306e-01 + u*( 2.8633188995618279e-02 + u*( -3.0850935191344057e-03 + u*( 2.6804897943923767e-04 + u*( -1.8646275300802761e-05 + u*( 1.0246499685157093e-06 + u*( -4.3510509586588465e-08 + u*( 1.3778914517450677e-09 + u*( -3.0648398393934102e-11 + u*( 4.2734150266051919e-13 + u*( -2.8114572543455210e-15))))))))))))))))));
                     case 19; w = 6.2041726508848370e+01 + u*( -1.1755274496413374e+02 + u*( 1.0517877181001440e+02 + u*( -5.9047731542464234e+01 + u*( 2.3308315082551665e+01 + u*( -6.8698191822257550e+00 + u*( 1.5668008661216632e+00 + u*( -2.8273098336030006e-01 + u*( 4.0921589696885537e-02 + u*( -4.7861508417409975e-03 + u*( 4.5342481658598943e-04 + u*( -3.4711947681224070e-05 + u*( 2.1314353839348120e-06 + u*( -1.0355151662841189e-07 + u*( 3.8929141589628539e-09 + u*( -1.0927478340948359e-10 + u*( 2.1567391462398081e-12 + u*( -2.6708843916282449e-14 + u*( 1.5619206968586228e-16))))))))))))))))));
                     otherwise; w = 0;
                end 
                case 20
                switch j 
                     case 1; w = 8.2206352466243317e+01 + u*( 1.5619206968586224e+02 + u*( 1.4057286271727605e+02 + u*( 7.9657955539789754e+01 + u*( 3.1863182215915902e+01 + u*( 9.5589546647747703e+00 + u*( 2.2304227551141134e+00 + u*( 4.1422136880690674e-01 + u*( 6.2133205321036014e-02 + u*( 7.5940584281266208e-03 + u*( 7.5940584281266221e-04 + u*( 6.2133205321036021e-05 + u*( 4.1422136880690667e-06 + u*( 2.2304227551141132e-07 + u*( 9.5589546647747719e-09 + u*( 3.1863182215915911e-10 + u*( 7.9657955539789751e-12 + u*( 1.4057286271727606e-13 + u*( 1.5619206968586228e-15 + u*( 8.2206352466243310e-18)))))))))))))))))));
                     case 2; w = -1.3989083239906907e+02 + u*( -3.1267976502979707e+02 + u*( -3.2829897199838337e+02 + u*( -2.1555764409599578e+02 + u*( -9.9343750955544365e+01 + u*( -3.4176689725711967e+01 + u*( -9.1084480127898537e+00 + u*( -1.9255456150462928e+00 + u*( -3.2782795865449738e-01 + u*( -4.5363630506822329e-02 + u*( -5.1247818166261115e-03 + u*( -4.7279294553703439e-04 + u*( -3.5481945634750973e-05 + u*( -2.1476339232897877e-06 + u*( -1.0333038813528232e-07 + u*( -3.8624549482133275e-09 + u*( -1.0817550362303450e-10 + u*( -2.1367075133025957e-12 + u*( -2.6552651846596591e-14 + u*( -1.5619206968586228e-16)))))))))))))))))));
                     case 3; w = 8.5205662588283204e+01 + u*( 2.2192441056516424e+02 + u*( 2.7313072554594856e+02 + u*( 2.1045505833123920e+02 + u*( 1.1366260025807316e+02 + u*( 4.5700691979394584e+01 + u*( 1.4189121651199548e+01 + u*( 3.4828187712369632e+00 + u*( 6.8624036377361364e-01 + u*( 1.0956347430858339e-01 + u*( 1.4241106285299615e-02 + u*( 1.5078092467053696e-03 + u*( 1.2956823705211605e-04 + u*( 8.9615129883262591e-06 + u*( 4.9180248212986317e-07 + u*( 2.0934747979501069e-08 + u*( 6.6673708786804028e-10 + u*( 1.4956952593118165e-11 + u*( 2.1085929407591400e-13 + u*( 1.4057286271727605e-15)))))))))))))))))));
                     case 4; w = -2.1619359279855296e+01 + u*( -6.8029220219783141e+01 + u*( -9.9666799748983934e+01 + u*( -9.1333414526563303e+01 + u*( -5.8787955660671024e+01 + u*( -2.8206689128638676e+01 + u*( -1.0446672051478188e+01 + u*( -3.0532081294734592e+00 + u*( -7.1433682923576258e-01 + u*( -1.3498174986765685e-01 + u*( -2.0693925739877578e-02 + u*( -2.5755061848088470e-03 + u*( -2.5931894690161883e-04 + u*( -2.0952885777345659e-05 + u*( -1.3396913198500503e-06 + u*( -6.6279242590970982e-08 + u*( -2.4480482896488189e-09 + u*( -6.3567048520752235e-11 + u*( -1.0355534220172666e-12 + u*( -7.9657955539789759e-15)))))))))))))))))));
                     case 5; w = 2.6508162588916582e+00 + u*( 8.8263356529156418e+00 + u*( 1.5616534060064156e+01 + u*( 1.7545289626426584e+01 + u*( 1.3797847107988865e+01 + u*( 8.0862122556912901e+00 + u*( 3.6672340424279035e+00 + u*( 1.3153818519736646e+00 + u*( 3.7781066612601855e-01 + u*( 8.7492739928261495e-02 + u*( 1.6385155892775494e-02 + u*( 2.4807322196438429e-03 + u*( 3.0248532025979103e-04 + u*( 2.9465445890986013e-05 + u*( 2.2616180850307837e-06 + u*( 1.3379350212463090e-07 + u*( 5.8883160735012602e-09 + u*( 1.8162013863072065e-10 + u*( 3.5049500437507498e-12 + u*( 3.1863182215915904e-14)))))))))))))))))));
                     case 6; w = 2.1984837742834515e-01 + u*( -4.1134229664484845e-01 + u*( -1.0112862491448782e+00 + u*( -1.2995733906769074e+00 + u*( -1.2780433056940144e+00 + u*( -9.5932199251840589e-01 + u*( -5.5401527340329981e-01 + u*( -2.5251075104935161e-01 + u*( -9.2557114780887226e-02 + u*( -2.7486050960093104e-02 + u*( -6.6106022848954398e-03 + u*( -1.2822100276113995e-03 + u*( -1.9924031270757459e-04 + u*( -2.4566545351653374e-05 + u*( -2.3696954500525920e-06 + u*( -1.7496073354759421e-07 + u*( -9.5493957101099954e-09 + u*( -3.6324027726144129e-10 + u*( -8.6030591982972938e-12 + u*( -9.5589546647747718e-14)))))))))))))))))));
                     case 7; w = 3.0743322578921345e-01 + u*( 4.6857330693378897e-03 + u*( -7.5223182288076368e-02 + u*( 2.6515954036820055e-02 + u*( 4.8046039019720627e-02 + u*( 3.5245016016857596e-02 + u*( 2.6148814908943802e-02 + u*( 1.6851147095615158e-02 + u*( 8.4535970234755793e-03 + u*( 3.3783332023510251e-03 + u*( 1.1054937557155955e-03 + u*( 2.9608234433176571e-04 + u*( 6.3808415949619674e-05 + u*( 1.0843860429122782e-05 + u*( 1.4242765978877094e-06 + u*( 1.4120360378076428e-07 + u*( 1.0210875372912409e-08 + u*( 5.0853638816601778e-10 + u*( 1.5612959285798793e-11 + u*( 2.2304227551141133e-13)))))))))))))))))));
                     case 8; w = 3.0669255881113333e-01 + u*( -5.1577921092119531e-06 + u*( -8.9295854872152886e-02 + u*( -6.5760844309566713e-05 + u*( 1.2603752511581481e-02 + u*( -1.9727049128118309e-04 + u*( -1.4174079307179292e-03 + u*( -2.1365751941307057e-04 + u*( -7.8805284038433034e-05 + u*( -9.7830700710216647e-05 + u*( -5.3227545304816835e-05 + u*( -1.9932555946527494e-05 + u*( -6.4171174455566264e-06 + u*( -1.7607224879601407e-06 + u*( -3.7637810455270773e-07 + u*( -5.8869140934837608e-08 + u*( -6.4618533533877470e-09 + u*( -4.7221236043987359e-10 + u*( -2.0711068440345339e-11 + u*( -4.1422136880690674e-13)))))))))))))))))));
                     case 9; w = 3.0669310173937742e-01 + u*( 2.6142622439794846e-11 + u*( -8.9272644690222294e-02 + u*( 1.3378815377167734e-09 + u*( 1.2735276875901773e-02 + u*( 1.6055212515791818e-08 + u*( -1.1872402931433713e-03 + u*( 6.9572620496017235e-08 + u*( 8.1490034987076083e-05 + u*( 1.2754980533441248e-07 + u*( -4.2484200470301013e-06 + u*( 1.0435893165697091e-07 + u*( 2.6185418050497516e-07 + u*( 3.7462180594899230e-08 + u*( 8.9471815662291461e-09 + u*( 5.3517400849852371e-09 + u*( 1.5657567740901071e-09 + u*( 2.3610618021993690e-10 + u*( 1.8639961596310804e-11 + u*( 6.2133205321036016e-13)))))))))))))))))));
                     case 10; w = 3.0669310173799164e-01 + u*( 2.2246205264930720e-14 + u*( -8.9272644926348912e-02 + u*( 2.0588008716879578e-14 + u*( 1.2735271524151709e-02 + u*( 1.8261191059231727e-15 + u*( -1.1872777553254197e-03 + u*( -1.0940685026790243e-16 + u*( 8.1385676055320891e-05 + u*( -7.0478683321661893e-18 + u*( -4.3759698523910332e-06 + u*( -6.8172394934309765e-20 + u*( 1.9228155940015619e-07 + u*( 2.2011711458253200e-22 + u*( -7.1080386887266927e-09 + u*( -1.2625388296862002e-23 + u*( 2.2782175284379820e-10 + u*( -7.6187687998305176e-26 + u*( -7.5940584281266223e-12 + u*( -7.5940584281266241e-13)))))))))))))))))));
                     case 11; w = 3.0669310173800551e-01 + u*( 3.9162968990165460e-15 + u*( -8.9272644926274736e-02 + u*( -3.4055178471838740e-15 + u*( 1.2735271524183097e-02 + u*( -2.4284276300226966e-15 + u*( -1.1872777553226567e-03 + u*( -1.1431414013200985e-16 + u*( 8.1385676055374992e-05 + u*( 2.8489170406823986e-18 + u*( -4.3759698523912500e-06 + u*( 2.7202748696062293e-20 + u*( 1.9228155940015847e-07 + u*( -9.8077499155761096e-22 + u*( -7.1080386887265298e-09 + u*( 1.2973674870568539e-23 + u*( 2.2782175284379833e-10 + u*( 7.6187687998305176e-26 + u*( -7.5940584281266223e-12 + u*( 7.5940584281266241e-13)))))))))))))))))));
                     case 12; w = 3.0669310173936576e-01 + u*( -2.6261939822794025e-11 + u*( -8.9272644690209860e-02 + u*( -1.3379686634301924e-09 + u*( 1.2735276875922579e-02 + u*( -1.6055226474877712e-08 + u*( -1.1872402931405162e-03 + u*( -6.9572621247729369e-08 + u*( 8.1490034987151760e-05 + u*( -1.2754980534682869e-07 + u*( -4.2484200470308416e-06 + u*( -1.0435893165706568e-07 + u*( 2.6185418050495186e-07 + u*( -3.7462180594900759e-08 + u*( 8.9471815662290303e-09 + u*( -5.3517400849852719e-09 + u*( 1.5657567740901059e-09 + u*( -2.3610618021993695e-10 + u*( 1.8639961596310804e-11 + u*( -6.2133205321036016e-13)))))))))))))))))));
                     case 13; w = 3.0669255881113278e-01 + u*( 5.1577919364389549e-06 + u*( -8.9295854872149402e-02 + u*( 6.5760844168664273e-05 + u*( 1.2603752511592769e-02 + u*( 1.9727049125895727e-04 + u*( -1.4174079307150936e-03 + u*( 2.1365751941192931e-04 + u*( -7.8805284038310709e-05 + u*( 9.7830700710192848e-05 + u*( -5.3227545304815629e-05 + u*( 1.9932555946527382e-05 + u*( -6.4171174455566197e-06 + u*( 1.7607224879601431e-06 + u*( -3.7637810455270746e-07 + u*( 5.8869140934837608e-08 + u*( -6.4618533533877428e-09 + u*( 4.7221236043987359e-10 + u*( -2.0711068440345335e-11 + u*( 4.1422136880690674e-13)))))))))))))))))));
                     case 14; w = 3.0743322578922561e-01 + u*( -4.6857330693277051e-03 + u*( -7.5223182288080739e-02 + u*( -2.6515954036843748e-02 + u*( 4.8046039019712126e-02 + u*( -3.5245016016861447e-02 + u*( 2.6148814908945262e-02 + u*( -1.6851147095614638e-02 + u*( 8.4535970234757684e-03 + u*( -3.3783332023509965e-03 + u*( 1.1054937557155972e-03 + u*( -2.9608234433176528e-04 + u*( 6.3808415949619647e-05 + u*( -1.0843860429122771e-05 + u*( 1.4242765978877087e-06 + u*( -1.4120360378076428e-07 + u*( 1.0210875372912405e-08 + u*( -5.0853638816601778e-10 + u*( 1.5612959285798793e-11 + u*( -2.2304227551141133e-13)))))))))))))))))));
                     case 15; w = 2.1984837742839228e-01 + u*( 4.1134229664491373e-01 + u*( -1.0112862491447703e+00 + u*( 1.2995733906769502e+00 + u*( -1.2780433056939822e+00 + u*( 9.5932199251841133e-01 + u*( -5.5401527340329693e-01 + u*( 2.5251075104935194e-01 + u*( -9.2557114780887087e-02 + u*( 2.7486050960093111e-02 + u*( -6.6106022848954425e-03 + u*( 1.2822100276113991e-03 + u*( -1.9924031270757462e-04 + u*( 2.4566545351653367e-05 + u*( -2.3696954500525920e-06 + u*( 1.7496073354759421e-07 + u*( -9.5493957101099938e-09 + u*( 3.6324027726144134e-10 + u*( -8.6030591982972938e-12 + u*( 9.5589546647747718e-14)))))))))))))))))));
                     case 16; w = 2.6508162588916697e+00 + u*( -8.8263356529157075e+00 + u*( 1.5616534060064199e+01 + u*( -1.7545289626426612e+01 + u*( 1.3797847107988895e+01 + u*( -8.0862122556912883e+00 + u*( 3.6672340424279080e+00 + u*( -1.3153818519736646e+00 + u*( 3.7781066612601849e-01 + u*( -8.7492739928261495e-02 + u*( 1.6385155892775491e-02 + u*( -2.4807322196438429e-03 + u*( 3.0248532025979103e-04 + u*( -2.9465445890986006e-05 + u*( 2.2616180850307832e-06 + u*( -1.3379350212463085e-07 + u*( 5.8883160735012585e-09 + u*( -1.8162013863072067e-10 + u*( 3.5049500437507498e-12 + u*( -3.1863182215915904e-14)))))))))))))))))));
                     case 17; w = -2.1619359279855235e+01 + u*( 6.8029220219783042e+01 + u*( -9.9666799748983721e+01 + u*( 9.1333414526563175e+01 + u*( -5.8787955660670967e+01 + u*( 2.8206689128638651e+01 + u*( -1.0446672051478179e+01 + u*( 3.0532081294734574e+00 + u*( -7.1433682923576247e-01 + u*( 1.3498174986765685e-01 + u*( -2.0693925739877567e-02 + u*( 2.5755061848088466e-03 + u*( -2.5931894690161877e-04 + u*( 2.0952885777345656e-05 + u*( -1.3396913198500501e-06 + u*( 6.6279242590970982e-08 + u*( -2.4480482896488189e-09 + u*( 6.3567048520752235e-11 + u*( -1.0355534220172666e-12 + u*( 7.9657955539789759e-15)))))))))))))))))));
                     case 18; w = 8.5205662588283261e+01 + u*( -2.2192441056516438e+02 + u*( 2.7313072554594862e+02 + u*( -2.1045505833123934e+02 + u*( 1.1366260025807316e+02 + u*( -4.5700691979394605e+01 + u*( 1.4189121651199553e+01 + u*( -3.4828187712369627e+00 + u*( 6.8624036377361342e-01 + u*( -1.0956347430858342e-01 + u*( 1.4241106285299615e-02 + u*( -1.5078092467053696e-03 + u*( 1.2956823705211602e-04 + u*( -8.9615129883262591e-06 + u*( 4.9180248212986317e-07 + u*( -2.0934747979501066e-08 + u*( 6.6673708786804028e-10 + u*( -1.4956952593118169e-11 + u*( 2.1085929407591405e-13 + u*( -1.4057286271727605e-15)))))))))))))))))));
                     case 19; w = -1.3989083239906907e+02 + u*( 3.1267976502979707e+02 + u*( -3.2829897199838337e+02 + u*( 2.1555764409599581e+02 + u*( -9.9343750955544365e+01 + u*( 3.4176689725711967e+01 + u*( -9.1084480127898519e+00 + u*( 1.9255456150462928e+00 + u*( -3.2782795865449732e-01 + u*( 4.5363630506822329e-02 + u*( -5.1247818166261115e-03 + u*( 4.7279294553703439e-04 + u*( -3.5481945634750979e-05 + u*( 2.1476339232897881e-06 + u*( -1.0333038813528232e-07 + u*( 3.8624549482133266e-09 + u*( -1.0817550362303447e-10 + u*( 2.1367075133025953e-12 + u*( -2.6552651846596591e-14 + u*( 1.5619206968586228e-16)))))))))))))))))));
                     case 20; w = 8.2206352466243317e+01 + u*( -1.5619206968586224e+02 + u*( 1.4057286271727605e+02 + u*( -7.9657955539789754e+01 + u*( 3.1863182215915902e+01 + u*( -9.5589546647747703e+00 + u*( 2.2304227551141134e+00 + u*( -4.1422136880690674e-01 + u*( 6.2133205321036014e-02 + u*( -7.5940584281266208e-03 + u*( 7.5940584281266221e-04 + u*( -6.2133205321036021e-05 + u*( 4.1422136880690667e-06 + u*( -2.2304227551141132e-07 + u*( 9.5589546647747719e-09 + u*( -3.1863182215915911e-10 + u*( 7.9657955539789751e-12 + u*( -1.4057286271727606e-13 + u*( 1.5619206968586228e-15 + u*( -8.2206352466243310e-18)))))))))))))))))));
                     otherwise; w = 0;
                end 
                end

                z(l)=w;
            end
        end%function
        
        function [M1_c,M1_r,M2_c,M2_r,M3_c,M3_r,M4_c,M4_r] = ToeplitzInv(c,r)
        % Computation of the inverse Toeplitz matrix of T(c,r) by means of the
        % Gohberg-Semencul formula
        % Call from outside class: infft.ToeplitzInv(c,r)
            if length(c)~=length(r)
                error('Computation only possible for quadratic Toeplitz matrix.')
            end

            n = size(c,1);
            e1 = zeros(n,1);
            en = zeros(n,1);
            e1(1) = 1;
            en(n) = 1;

            % Prevent warning for tiny difference between r(1) and c(1)
            r(1) = c(1);

            % Compute the relevant vectors
            x = toeplitz(c,r)\e1;
            y = toeplitz(c,r)\en;

            % Compute the components of the Gohberg-Semencul formula
            M1_c = x; % First column of M1
            M1_r = [x(1),zeros(1,n-1)]; % First row of M1
            M2_c = [y(n);zeros(n-1,1)]; % First column of M2
            M2_r = flip(y).'; % First row of M2
            M3_c = [0;y(1:n-1)]; % First column of M2
            M3_r = zeros(1,n); % First row of M3
            M4_c = zeros(n,1); % First column of M3
            M4_r = [0,flip(x(2:n)).']; % First row of M4
        end%function
        
        function D = ToeplitzDiag(c,r)
        % Diagonalization of circular embedding of Toeplitz matrix T(c,r) via FFT
        % Call from outside class: infft.ToeplitzDiag(c,r)
            n = length(c);
            C = [c;0;r(n:-1:2).'];
            D = fft(C);
            D = D/sqrt(2*n);
        end%function
        
        function u = ToeplitzMul(D,v)
        % Multiplication of Toeplitz matrix T(c,r) with vector v 
        % via given diagonalized embedding
        % Call from outside class: infft.ToeplitzMul(D,v)
            n = length(v);
            v = [v(:);zeros(n,1)]; % zero-padding

            % FFT
            w = fft(v);
            w = w/sqrt(2*n);

            % Multiplication with Diagonalization and iFFT
            u = ifft(D.*w);
            u = 2*n*u;

            % Discard last n elements
            u = u(1:n);
        end%function
        
        function u = ToeplitzInvMul(x0,D1,D2,D3,D4,v)
        % Multiplication of inverse of Toeplitz matrix T(c,r) with vector v
        % via given diagonalized embedding
        % Call from outside class: infft.ToeplitzInvMul(x1,D1,D2,D3,D4,v)
            u1 = infft.ToeplitzMul(D2,v);
            u1 = infft.ToeplitzMul(D1,u1);
            u2 = infft.ToeplitzMul(D4,v);
            u2 = infft.ToeplitzMul(D3,u2);

            % Apply Gohberg-Semencul formula
            u = (u1-u2)/x0;
        end%function
    end %methods
end %classdef