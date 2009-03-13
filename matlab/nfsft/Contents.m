%
% Files
%   f_hat_index               - Return indices of Fourier coefficients
%   ndsft_adjoint             - Adjoint discrete spherical Fourier transform (direct alg.)
%   ndsft_trafo               - Discrete spherical Fourier transform (direct algorithm)
%   nfsft_adjoint             - Adjoint fast spherical Fourier transform
%   nfsft_finalize            - Finalize plan
%   nfsft_forget              - Forget precomputed data
%   nfsft_get_f               - Get function values from plan
%   nfsft_get_f_hat           - Get Fourier coefficients in a matrix from plan
%   nfsft_get_x               - Get nodes from plan
%   nfsft_init                - Plan initialization
%   nfsft_init_advanced       - Advanced plan initialization routine
%   nfsft_init_guru           - Guru plan initialization
%   nfsft_precompute          - Node-independent precomputation (for FPT)
%   nfsft_precompute_x        - Node-dependent precomputation (for NFFT)
%   nfsft_set_f               - Set function values in plan
%   nfsft_set_f_hat           - Set Fourier coefficients in plan from a matrix
%   nfsft_set_x               - Set nodes in plan
%   nfsft_trafo               - Fast spherical Fourier transform
%   NFSFT_NORMALIZED          - Flag for using L^2-normalized Spherical harmonics
%   NFSFT_NO_DIRECT_ALGORITHM - Flag for not precomputing for the direct algorithm
%   NFSFT_NO_FAST_ALGORITHM   - Flag for not precomputing for the fast algorithm
%   NFSFT_PRESERVE_F_HAT      - Flag for NFSFT not destroying input in f_hat
%   NFSFT_USE_DPT             - Flag for using the DPT algorithm internally
%   NFSFT_USE_NDFT            - Flag for using the NDFT algorithm internally
%   cc                        - Clenshaw-Curtis interpolatory quadrature rule
%   gl                        - Gauss-Legendre interpolatory quadrature rule
%   nfsft_get_f_hat_linear    - Get Fourier coefficients in a linear vector from plan
%   nfsft_set_f_hat_linear    - Set Fourier coefficients in plan from a linear vector
%   projection                - Example program: Project a function onto spherical harmonics.
%   simple_test               - Example program: Basic usage principles
%   nfsftmex                  - Gateway routine to the NFSFT module
