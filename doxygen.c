/** \mainpage
 *
 * \section section_intro_sec Introduction
 *
 * Fast Fourier transforms (FFTs) belong to the '10 algorithms with the
 * greatest influence on the development and practice of science and
 * engineering in the 20th century' (J. Dongarra, S. Sullivan).
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
 * We believe that our library, which is free software and based on  FFTW3,
 * should become the NFFT library of choice for most applications.
 */
