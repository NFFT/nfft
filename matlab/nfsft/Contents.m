%
% Files
%   f_hat_index         - Return indices of Fourier coefficients
%   ndsft_adjoint       - Adjoint fast spherical Fourier transformation (direct algorithm)
%   ndsft_trafo         - Fast spherical Fourier transformation (direct algorithm)
%   nfsft               - Gateway function to NFSFT module from NFFT3
%   nfsft_adjoint       - Adjoint fast spherical Fourier transformation
%   nfsft_finalize      - Finalize plan
%   nfsft_forget        - Forget precomputed data
%   nfsft_get_f         - Get function values from plan
%   nfsft_get_f_hat     - Get Fourier coefficients from plan
%   nfsft_get_x         - Get nodes from plan
%   nfsft_init          - Initialise plans
%   nfsft_init_advanced - Initialise plans (advanced)
%   nfsft_init_guru     - Initialise plans (guru)
%   nfsft_precompute    - Precompute
%   nfsft_precompute_x  - Precompute for x
%   nfsft_set_f         - Set function values in plan
%   nfsft_set_f_hat     - Set Fourier coefficients in plan
%   nfsft_set_x         - Set nodes in plan
%   nfsft_trafo         - fast spherical Fourier transformation
%   NFSFT_NORMALIZED    - Flag for using L^2-normalized Spherical harmonics
%   NFSFT_NO_DIRECT_ALGORITHM - Flag for not precomputing for the direct algorithm
%   NFSFT_NO_FAST_ALGORITHM   - Flag for not precomputing for the fast algorithm
%   NFSFT_PRESERVE_F_HAT      - Flags for NFSFT not destroying input in f_hat
%   NFSFT_USE_DPT             - Flag for using the DPT algorithm internally.
%   NFSFT_USE_NDFT            - Flag for using the NDFT algorithm internally.
