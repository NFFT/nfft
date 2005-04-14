#ifndef DIRECT_H
#define DIRECT_H

#include "api.h"

/**
 * Direct spherical Fourier transform
 */
void direct_trafo(double *angles, int D, complex **f_hat, complex *f, 
  int M, int N, struct nfsft_transform_wisdom *tw, struct nfsft_wisdom *wisdom);
void direct_trafo_adjoint(double *angles, int D, complex **f_hat, complex *f, 
  int M, int N, struct nfsft_transform_wisdom *tw, struct nfsft_wisdom *wisdom);
#endif

