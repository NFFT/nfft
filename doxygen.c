/** \mainpage
 *
 * \section section_intro_sec Introduction
 *
 * Fast Fourier transforms (FFTs) belong to the '10 algorithms with the
 * greatest influence on the development and practice of science and
 * engineering in the 20th century'.
 * The classic algorithm computes the discrete Fourier transform
 * \f[
 *   f_j= \sum_{k=-\frac{N}{2}}^{\frac{N}{2}-1} \hat{f}_{k}
 * {\rm e}^{2\pi{\rm i}\frac{kj}{N}}
 * \f]
 * for \f$j=-\frac{N}{2},\hdots,\frac{N}{2}-1\f$ and given complex
 * coefficients \f$\hat{f}_{k}\in\mathbb{C}\f$.
 * Using a divide and conquer approach, the number of floating point
 * operations is reduced from \f${\cal O}(N^2)\f$ for a straightforward
 * computation to only \f${\cal O}(N\log N)\f$.
 * In conjunction with publicly available efficient implementations the fast
 * Fourier transform has become of great importance in scientific computing.

 * However, two shortcomings of traditional schemes are the need for
 * equispaced sampling and the restriction to the system of complex
 * exponential functions.
 * The NFFT is a C subroutine library for computing the nonequispaced discrete
 * Fourier transform (NDFT) and its generalisations in one or more dimensions,
 * of arbitrary input size, and of complex data.
 *
 * More precisely,we collect the possible frequencies 
 * \f$\mathbf{k}\in\mathbb{Z}^d\f$ in the multi-index set
 *\f[
 * I_{\mathbf{N}} := \left\{ \mathbf{k}=\left(k_t\right)_{t=0,\hdots,d-1} 
 *  \in \mathbb{Z}^d: -
 *   \frac{N_t}{2} \le k_t < \frac{N_t}{2} ,\;t=0,\hdots,d-1\right\},
 *\f]
 * where \f$\mathbf{N}=\left(N_t\right)_{t=0,\hdots,d-1}\f$ is the
 * multibandlimit, i.e., \f$N_t\in 2\mathbb{N}\f$.
 * For a finite number of given Fourier coefficients 
 * \f$\hat f_{\mathbf{k}} \in \mathbb{C}\f$, 
 * \f$\mathbf{k}\in I_{\mathbf{N}}\f$, we consider the 
 * fast evaluation of the trigonometric polynomial 
 * \f[
 *  f\left(\mathbf{x}\right) 
 *  := \sum_{ \mathbf{k}\in I_{ N}} \hat{f}_{\mathbf{ k}} 
 *  {\rm e}^{-2\pi{\rm i}\mathbf{k}\mathbf{ x}}
 * \f]
 * at given nonequispaced nodes \f$\mathbf{x}_j \in \mathbb{T}^d\f$,
 * \f$j=0,\ldots, M-1\f$, from the 
 * \f$ d\f$-dimensional torus as well as the 
 * adjoint problem, the fast evaluation of sums of the form
 * \f[
 *  \hat h_{\mathbf{k}} := \sum_{j=0}^{M-1} {f}_{j} 
 *  {\rm e}^{2\pi{\rm i}\mathbf{k}\mathbf{ x}_j}.
 * \f]
 *
 * The generalisations of the NFFT include
 *    - NNFFT - nonequispaced in time and frequency fast Fourier transform,
 *    - NFCT/NFST - nonequispaced fast (co)sine transform,
 *    - NSFFT - nonequispaced sparse fast Fourier transform,
 *    - FPT - fast polynomial transform,
 *    - NFSFT - nonequispaced fast spherical Fourier transform.
 *  
 *
 *
 * Furthermore, we consider the inversion of the above transforms by 
 * iterative methods.
 */
