% NFCT
%
% Files
%   FFT_OUT_OF_PLACE    - FFT flag
%   FFTW_ESTIMATE       - FFT flag
%   FFTW_MEASURE        - FFT flag
%   FG_PSI              - Precompuation flag
%   ndct_adjoint        - Adjoint nonequispaced discrete cosine transform (direct algorithm)
%   ndct_trafo          - Nonequispaced discrete cosine transform (direct algorithm)
%   nfct_adjoint        - Adjoint nonequispaced fast cosine transform
%   nfct_finalize       - Finalize plan
%   nfct_get_f          - Get function values from plan
%   nfct_get_f_hat      - Get Fourier coefficients from plan
%   nfct_get_x          - Get nodes from plan
%   nfct_init_1d        - Initialise plans
%   nfct_init_2d        - Initialise plans
%   nfct_init_3d        - Initialise plans
%   nfct_init_guru      - Initialise plans, no error handling
%   nfct_precompute_psi - Precompute psi, dependent on nodes x
%   nfct_set_f          - Set function values in plan
%   nfct_set_f_hat      - Set Fourier coefficients in plan
%   nfct_set_x          - Set nodes in plan
%   nfct_trafo          - nonequispaced fast cosine transform
%   nfctmex             - Gateway function to NFCT module from NFFT3
%   PRE_FG_PSI          - Precomputation flag
%   PRE_FULL_PSI        - Precomputation flag
%   PRE_LIN_PSI         - Precomputation flag
%   PRE_PHI_HUT         - Precomputation flag
%   PRE_PSI             - Precomputation flag
%   simple_test         - Example program: Basic usage principles
