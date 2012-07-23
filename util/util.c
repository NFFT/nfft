/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* $Id$ */

/** Sources for utilities.
 *  functions for vectors, window functions, ...
 *  (c) if not stated otherwise: Daniel Potts, Stefan Kunis
 */
#include "config.h"

#include "infft.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>
#include "cstripack.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "nfft3.h"
#include "nfft3util.h"
#include "infft.h"

double nfft_elapsed_seconds(ticks t1, ticks t0)
{
  UNUSED(t1);
  UNUSED(t0);
  return elapsed(t1,t0) / TICKS_PER_SECOND;
}

/** Computes integer /f$\prod_{t=0}^{d-1} v_t/f$.
 */
int nfft_prod_int(int *vec, int d)
{
  int t, prod;

  prod=1;
  for(t=0; t<d; t++)
    prod *= vec[t];

  return prod;
}

/** Computes integer /f$\prod_{t=0}^{d-1} v_t-a/f$.
 */
int nfst_prod_minus_a_int(int *vec, int a, int d)
{
  int t, prod;

  prod=1;
  for(t=0; t<d; t++)
    prod *= vec[t]-a;

  return prod;
}

/** Computes /f$\sum_{t=0}^{d-1} i_t \prod_{t'=t+1}^{d-1} N_{t'}/f$.
 */
int nfft_plain_loop(int *idx,int *N,int d)
{
  int t,sum;

  sum = idx[0];
  for (t = 1; t < d; t++)
    sum = sum * N[t] + idx[t];

  return sum;
}

/** Computes double /f$\prod_{t=0}^{d-1} v_t/f$.
 */
double nfft_prod_real(double *vec,int d)
{
  int t;
  double prod;

  prod=1.0;
  for(t=0; t<d; t++)
    prod*=vec[t];

  return prod;
}

/** Computes the inner/dot product \f$x^H x\f$.
 */
double nfft_dot_complex(double _Complex *x, int n)
{
  int k;
  double dot;

  for(k=0,dot=0; k<n; k++)
    dot+=conj(x[k])*x[k];

  return dot;
}

/** Computes the inner/dot product \f$x^H x\f$.
 */
double nfft_dot_double(double *x, int n)
{
  int k;
  double dot;

  for(k=0,dot=0; k<n; k++)
    dot+=x[k]*x[k];

  return dot;
}


/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double nfft_dot_w_complex(double _Complex *x, double *w, int n)
{
  int k;
  double dot;

  for(k=0,dot=0.0; k<n; k++)
    dot+=w[k]*conj(x[k])*x[k];

  return dot;
}

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double nfft_dot_w_double(double *x, double *w, int n)
{
  int k;
  double dot;

  for(k=0,dot=0.0; k<n; k++)
    dot+=w[k]*x[k]*x[k];

  return dot;
}


/** Computes the weighted inner/dot product
    \f$x^H (w\odot w2\odot w2 \odot x)\f$.
 */
double nfft_dot_w_w2_complex(double _Complex *x, double *w, double *w2, int n)
{
  int k;
  double dot;

  for(k=0,dot=0.0; k<n; k++)
    dot+=w[k]*w2[k]*w2[k]*conj(x[k])*x[k];

  return dot;
}

/** Computes the weighted inner/dot product
    \f$x^H (w2\odot w2 \odot x)\f$.
 */
double nfft_dot_w2_complex(double _Complex *x, double *w2, int n)
{
  int k;
  double dot;

  for(k=0,dot=0.0; k<n; k++)
    dot+=w2[k]*w2[k]*conj(x[k])*x[k];

  return dot;
}

/** Copies \f$x \leftarrow y\f$.
 */
void nfft_cp_complex(double _Complex *x, double _Complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=y[k];
}

/** Copies \f$x \leftarrow y\f$.
 */
void nfft_cp_double(double *x, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=y[k];
}

/** Copies \f$x \leftarrow a y\f$.
 */
void nfft_cp_a_complex(double _Complex *x, double a, double _Complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*y[k];
}

/** Copies \f$x \leftarrow a y\f$.
 */
void nfft_cp_a_double(double *x, double a, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*y[k];
}


/** Copies \f$x \leftarrow w\odot y\f$.
 */
void nfft_cp_w_complex(double _Complex *x, double *w, double _Complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=w[k]*y[k];
}

/** Copies \f$x \leftarrow w\odot y\f$.
 */
void nfft_cp_w_double(double *x, double *w, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=w[k]*y[k];
}



/** Updates \f$x \leftarrow a x + y\f$.
 */
void nfft_upd_axpy_complex(double _Complex *x, double a, double _Complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+y[k];
}

/** Updates \f$x \leftarrow a x + y\f$.
 */
void nfft_upd_axpy_double(double *x, double a, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+y[k];
}


/** Updates \f$x \leftarrow x + a y\f$.
 */
void nfft_upd_xpay_complex(double _Complex *x, double a, double _Complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]+=a*y[k];
}

/** Updates \f$x \leftarrow x + a y\f$.
 */
void nfft_upd_xpay_double(double *x, double a, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]+=a*y[k];
}



/** Updates \f$x \leftarrow a x + b y\f$.
 */
void nfft_upd_axpby_complex(double _Complex *x, double a, double _Complex *y, double b, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+b*y[k];
}

/** Updates \f$x \leftarrow a x + b y\f$.
 */
void nfft_upd_axpby_double(double *x, double a, double *y, double b, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+b*y[k];
}


/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void nfft_upd_xpawy_complex(double _Complex *x, double a, double *w, double _Complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]+=a*w[k]*y[k];
}

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void nfft_upd_xpawy_double(double *x, double a, double *w, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]+=a*w[k]*y[k];
}



/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void nfft_upd_axpwy_complex(double _Complex *x, double a, double *w, double _Complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+w[k]*y[k];
}

/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void nfft_upd_axpwy_double(double *x, double a, double *w, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+w[k]*y[k];
}


void nfft_fftshift_complex(double _Complex *x, int d, int* N)
{
  int d_pre, d_act, d_post;
  int N_pre, N_act, N_post;
  int k_pre, k_act, k_post;
  int k,k_swap;

  double _Complex x_swap;

  for(d_act=0;d_act<d;d_act++)
    {
      for(d_pre=0, N_pre=1;d_pre<d_act;d_pre++)
  N_pre*=N[d_pre];

      N_act=N[d_act];

      for(d_post=d_act+1, N_post=1;d_post<d;d_post++)
  N_post*=N[d_post];

      for(k_pre=0;k_pre<N_pre;k_pre++)
  for(k_act=0;k_act<N_act/2;k_act++)
    for(k_post=0;k_post<N_post;k_post++)
      {
        k=(k_pre*N_act+k_act)*N_post+k_post;
        k_swap=(k_pre*N_act+k_act+N_act/2)*N_post+k_post;

        x_swap=x[k];
        x[k]=x[k_swap];
        x[k_swap]=x_swap;
      }
    }
}

/** vector print
 */
void nfft_vpr_int(int *x, int n, char *text)
{
  int k;

  if(text!=NULL)
  {
      printf ("\n %s, adr=%p\n", text, (void*)x);
      for (k=0; k<n; k++)
      {
    if (k%8==0)
        printf("%6d.\t", k);
    printf("%d,", x[k]);
    if (k%8==7)
        printf("\n");
      }
      if (n%8!=0)
        printf("\n");
  }
  else
      for (k=0; k<n; k++)
    printf("%d,\n", x[k]);
  fflush(stdout);
}

/** Print real vector to standard output. */
void X(vpr_double)(R *x, const int n, const char *text)
{
  int k;

  if (x == NULL)
  {
    printf("null pointer\n");
    fflush(stdout);
    exit(-1);
  }

  if (text != NULL)
  {
    printf ("\n %s, adr=%p\n", text, (void*)x);

    for (k = 0; k < n; k++)
    {
      if (k%8 == 0)
        printf("%6d.\t", k);

      printf("%+.1" FE ",", x[k]);

      if (k%8 == 7)
        printf("\n");
    }

    if (n%8 != 0)
      printf("\n");
  }
  else
    for (k = 0; k < n; k++)
      printf("%+" FE ",\n", x[k]);

  fflush(stdout);
}

/** Print complex vector to standard output. */
void X(vpr_complex)(C *x, const int n, const char *text)
{
  int k;

  if(text != NULL)
  {
    printf("\n %s, adr=%p\n", text, (void*)x);
    for (k = 0; k < n; k++)
    {
      if (k%4 == 0)
        printf("%6d.\t", k);

      printf("%+.1" FE "%+.1" FE "i,", CREAL(x[k]), CIMAG(x[k]));

      if (k%4==3)
        printf("\n");
    }
    if (n%4!=0)
      printf("\n");
  }
  else
    for (k = 0; k < n; k++)
      printf("%+" FE "%+" FE "i,\n", CREAL(x[k]), CIMAG(x[k]));

  fflush(stdout);
}

/** Compute non periodic voronoi weights for ordered nodes x_j */
void X(voronoi_weights_1d)(R *w, R *x, const int M)
{
  int j;

  w[0] = (x[1]-x[0])/K(2.0);

  for(j = 1; j < M-1; j++)
    w[j] = (x[j+1]-x[j-1])/K(2.0);

  w[M-1] = (x[M-1]-x[M-2])/K(2.0);
}

void nfft_voronoi_weights_S2(double *w, double *xi, int M)
{
  double *x;
  double *y;
  double *z;
  int j;
  int k;
  int el;
  int Mlocal = M;
  int lnew;
  int ier;
  int *list;
  int *lptr;
  int *lend;
  int *near;
  int *next;
  double  *dist;
  int *ltri;
  int *listc;
  int nb;
  double *xc;
  double *yc;
  double *zc;
  double *rc;
  double *vr;
  int lp;
  int lpl;
  int kv;
  double a;

  /* Allocate memory for auxilliary arrays. */
  x = (double*)nfft_malloc(M * sizeof(double));
  y = (double*)nfft_malloc(M * sizeof(double));
  z = (double*)nfft_malloc(M * sizeof(double));

  list = (int*)nfft_malloc((6*M-12+1)*sizeof(int));
  lptr = (int*)nfft_malloc((6*M-12+1)*sizeof(int));
  lend = (int*)nfft_malloc((M+1)*sizeof(int));
  near = (int*)nfft_malloc((M+1)*sizeof(int));
  next = (int*)nfft_malloc((M+1)*sizeof(int));
  dist = (double*)nfft_malloc((M+1)*sizeof(double));
  ltri = (int*)nfft_malloc((6*M+1)*sizeof(int));
  listc = (int*)nfft_malloc((6*M-12+1)*sizeof(int));
  xc = (double*)nfft_malloc((2*M-4+1)*sizeof(double));
  yc = (double*)nfft_malloc((2*M-4+1)*sizeof(double));
  zc = (double*)nfft_malloc((2*M-4+1)*sizeof(double));
  rc = (double*)nfft_malloc((2*M-4+1)*sizeof(double));
  vr = (double*)nfft_malloc(3*(2*M-4+1)*sizeof(double));

  /* Convert from spherical Coordinates in [0,1/2]x[-1/2,1/2) to Cartesian
   * coordinates. */
  for (k = 0; k < M; k++)
  {
    x[k] = sin(2.0*PI*xi[2*k+1])*cos(2.0*PI*xi[2*k]);
    y[k] = sin(2.0*PI*xi[2*k+1])*sin(2.0*PI*xi[2*k]);
    z[k] = cos(2.0*PI*xi[2*k+1]);
  }

  /* Generate Delaunay triangulation. */
  trmesh_(&Mlocal, x, y, z, list, lptr, lend, &lnew, near, next, dist, &ier);

  /* Check error flag. */
  if (ier == 0)
  {
    /* Get Voronoi vertices. */
    crlist_(&Mlocal, &Mlocal, x, y, z, list, lend, lptr, &lnew, ltri, listc, &nb, xc,
      yc, zc, rc, &ier);

    if (ier == 0)
    {
      /* Calcuate sizes of Voronoi regions. */
      for (k = 0; k < M; k++)
      {
        /* Get last neighbour index. */
        lpl = lend[k];
        lp = lpl;

        j = 0;
        vr[3*j] = x[k];
        vr[3*j+1] = y[k];
        vr[3*j+2] = z[k];

        do
        {
          j++;
          /* Get next neighbour. */
          lp = lptr[lp-1];
          kv = listc[lp-1];
          vr[3*j] = xc[kv-1];
          vr[3*j+1] = yc[kv-1];
          vr[3*j+2] = zc[kv-1];
          /* fprintf(stderr, "lp = %ld\t", lp); */
        } while (lp != lpl);

        a = 0;
        for (el = 0; el < j; el++)
        {
          a += areas_(vr, &vr[3*(el+1)],&vr[3*(((el+1)%j)+1)]);
        }

        w[k] = a;
      }
    }
  }

  /* Deallocate memory. */
  nfft_free(x);
  nfft_free(y);
  nfft_free(z);

  nfft_free(list);
  nfft_free(lptr);
  nfft_free(lend);
  nfft_free(near);
  nfft_free(next);
  nfft_free(dist);
  nfft_free(ltri);
  nfft_free(listc);
  nfft_free(xc);
  nfft_free(yc);
  nfft_free(zc);
  nfft_free(rc);
  nfft_free(vr);
}

/**
 * Compute damping factor for modified Fejer kernel:
 * /f$\frac{2}{N}\left(1-\frac{\left|2k+1\right|}{N}\right)/f$
 */
R X(modified_fejer)(const int N, const int kk)
{
  return (K(2.0)/((R)(N*N))*(K(1.0)-FABS(K(2.0)*kk+K(1.0))/((R)N)));
}

/** Compute damping factor for modified Jackson kernel. */
R X(modified_jackson2)(const int N, const int kk)
{
  int kj;
  const R n=(N/K(2.0)+K(1.0))/K(2.0);
  R result, k;

  for (result = K(0.0), kj = kk; kj <= kk+1; kj++)
  {
    k = ABS(kj);

    if(k/n < K(1.0))
      result += K(1.0) - (K(3.0)*k + K(6.0)*n*POW(k,K(2.0))
        - K(3.0)*POW(k,K(3.0)))/(K(2.0)*n*(K(2.0)*POW(n,K(2.0))+K(1.0)));
    else
      result+= (K(2.0)*n-k)*(POW(2*n-k,K(2.0))-K(1.0))/(K(2.0)
        *n*(K(2.0)*POW(n,K(2.0))+K(1.0)));
  }

  return result;
}

/** Compute damping factor for modified generalised Jackson kernel. */
R X(modified_jackson4)(const int N, const int kk)
{
  int kj;
  const R n = (N/K(2.0)+K(3.0))/K(4.0), normalisation = (K(2416.0)*POW(n,K(7.0))
    + K(1120.0)*POW(n,K(5.0)) + K(784.0)*POW(n,K(3.0)) + K(720.0)*n);
  R result, k;

  for (result = K(0.0), kj = kk; kj <= kk + 1; kj++)
  {
    k = ABS(kj);

    if (k/n < K(1.0))
      result += K(1.0) - (K(1260.0)*k + (K(1680.0)*POW(n, K(5.0))
        + K(2240.0)*POW(n, K(3.0)) + K(2940.0)*n)*POW(k, K(2.0))
        - K(1715.0)*POW(k, K(3.0)) - (K(560.0)*POW(n, K(3.0))
        + K(1400.0)*n)*POW(k, K(4.0)) + K(490.0)*POW(k, K(5.0))
        + K(140.0)*n*POW(k, K(6.0)) - K(35.0)*POW(k,K(7.0)))/normalisation;

    if ((K(1.0) <= k/n) && (k/n < K(2.0)))
      result += ((K(2472.0)*POW(n, K(7.0)) + K(336.0)*POW(n, K(5.0))
        + K(3528.0)*POW(n, K(3.0)) - K(1296.0)*n) - (K(392.0)*POW(n, K(6.0))
        - K(3920.0)*POW(n, K(4.0)) + K(8232.0)*POW(n, K(2.0)) - K(756.0))*k
        - (K(504.0)*POW(n, K(5.0)) + K(10080.0)*POW(n, K(3.0))
        - K(5292.0)*n)*POW(k, K(2.0)) - (K(1960.0)*POW(n, K(4.0))
        - K(7840.0)*POW(n, K(2.0)) + K(1029.0))*POW(k, K(3.0))
        + (K(2520.0)*POW(n, K(3.0)) - K(2520.0)*n) * POW(k, K(4.0))
        - (K(1176.0)*POW(n, K(2.0)) - K(294.0)) * POW(k, K(5.0))
        + K(252.0)*n*POW(k, K(6.0)) - K(21.0)*POW(k, K(7.0)))/normalisation;

    if ((K(2.0) <= k/n) && (k/n < K(3.0)))
      result += (-(K(1112.0)*POW(n, K(7.0)) - K(12880.0)*POW(n, K(5.0))
        + K(7448.0)*POW(n, K(3.0)) - K(720.0)*n) + (K(12152.0)*POW(n, K(6.0))
        - K(27440.0)*POW(n, K(4.0)) + K(8232.0)*POW(n, K(2.0)) - K(252.0))*k
        - (K(19320.0)*POW(n, K(5.0)) - K(21280.0)*POW(n, K(3.0))
        + K(2940.0)*n)*POW(k, K(2.0)) + (K(13720.0)*POW(n, K(4.0))
        - K(7840.0)*POW(n, K(2.0)) + K(343.0))*POW(k, K(3.0))
        - (K(5320.0)*POW(n, K(3.0)) - K(1400.0)*n)*POW(k, K(4.0))
        + (K(1176.0)*POW(n, K(2.0)) - K(98.0))*POW(k, K(5.0))
        - K(140.0)*n*POW(k, K(6.0)) + K(7.0) * POW(k, K(7.0)))/normalisation;

    if ((K(3.0) <= k/n) && (k/n < K(4.0)))
      result += ((4*n-k)*(POW(4*n-k, K(2.0)) - K(1.0))*(POW(4*n-k, K(2.0))
        - K(4.0))*(POW(4*n-k, K(2.0)) - K(9.0)))/normalisation;
  }

  return result;
}

/** Compute damping factor for modified Sobolev kernel. */
R X(modified_sobolev)(const R mu, const int kk)
{
  R result;
  int kj, k;

  for (result = K(0.0), kj = kk; kj <= kk+1; kj++)
  {
    k = ABS(kj);
    if (k == 0)
      result += K(1.0);
    else
      result += POW((double)k,-K(2.0)*mu);
  }

  return result;
}

/** Comput damping factor for modified multiquadric kernel. */
R X(modified_multiquadric)(const R mu, const R c, const int kk)
{
  R result;
  int kj, k;

  for (result = K(0.0), kj = kk; kj <= kk+1; kj++)
    {
      k = ABS(kj);
      result += POW((double)(k*k + c*c), -mu);
    }

  return result;
}

static inline int scaled_modified_bessel_i_series(const R x, const R alpha,
  const int nb, const int ize, R *b)
{
  const R enmten = K(4.0)*nfft_float_property(NFFT_R_MIN);
  R tempa = K(1.0), empal = K(1.0) + alpha, halfx = K(0.0), tempb = K(0.0);
  int n, ncalc = nb;

  if (enmten < x)
    halfx = x/K(2.0);

  if (alpha != K(0.0))
    tempa = POW(halfx, alpha)/TGAMMA(empal);

  if (ize == 2)
    tempa *= EXP(-x);

  if (K(1.0) < x + K(1.0))
    tempb = halfx*halfx;

  b[0] = tempa + tempa*tempb/empal;

  if (x != K(0.0) && b[0] == K(0.0))
    ncalc = 0;

  if (nb == 1)
    return ncalc;

  if (K(0.0) < x)
  {
    R tempc = halfx, tover = (enmten + enmten)/x;

    if (tempb != K(0.0))
      tover = enmten/tempb;

    for (n = 1; n < nb; n++)
    {
      tempa /= empal;
      empal += K(1.0);
      tempa *= tempc;

      if (tempa <= tover*empal)
        tempa = K(0.0);

      b[n] = tempa + tempa*tempb/empal;

      if (b[n] == K(0.0) && n < ncalc)
        ncalc = n;
    }
  }
  else
    for (n = 1; n < nb; n++)
      b[n] = K(0.0);

  return ncalc;
}

static inline void scaled_modified_bessel_i_normalize(const R x,
  const R alpha, const int nb, const int ize, R *b, const R sum_)
{
  const R enmten = K(4.0)*nfft_float_property(NFFT_R_MIN);
  R sum = sum_, tempa;
  int n;

  /* Normalize, i.e., divide all b[n] by sum */
  if (alpha != K(0.0))
    sum = sum * TGAMMA(K(1.0) + alpha) * POW(x/K(2.0), -alpha);

  if (ize == 1)
    sum *= EXP(-x);

  tempa = enmten;

  if (K(1.0) < sum)
    tempa *= sum;

  for (n = 1; n <= nb; n++)
  {
    if (b[n-1] < tempa)
      b[n-1] = K(0.0);

    b[n-1] /= sum;
  }
}

/**
 * Calculates the modified bessel function \f$I_{n+\alpha}(x)\f$, possibly
 * scaled by \f$\mathrm{e}^{-x}\f$, for real non-negative \f$x,alpha\f$ with
 * \f$0 \le \alpha < 1\f$, and \f$n=0,1,\ldots,nb-1\f$.
 *
 * \arg[in] \c x non-negative real number in \f$I_{n+\alpha}(x)\f$
 * \arg[in] \c alpha non-negative real number with \f$0 \le \alpha < 1\f$ in
 *   \f$I_{n+\alpha}(x)\f$
 * \arg[in] \c nb number of functions to be calculated
 * \arg[in] \c ize switch between no scaling (\c ize = 1) and exponential
 *   scaling (\c ize = 2)
 * \arg[out] \c b real output vector to contain \f$I_{n+\alpha}(x)\f$,
 *   \f$n=0,1,\ldots,nb-1\f$
 * \return error indicator. Only if this value is identical to \c nb, then all
 *   values in \c b have been calculated to full accuracy. If not, errors are
 *   indicated using the following scheme:
 *   - ncalc < 0: At least one of the arguments was out of range (e.g. nb <= 0,
 *     ize neither equals 1 nor 2, \f$|x| \ge exparg\f$). In this case, the
 *     output vector b is not calculated and \c ncalc is set to
 *     \f$\min(nb,0)-1\f$.
 *   - 0 < ncalc < nb: Not all requested functions could be calculated to full
 *     accuracy. This can occur when nb is much larger than |x|. in this case,
 *     the values \f$I_{n+\alpha}(x)\f$ are calculated to full accuracy for
 *     \f$n=0,1,\ldots,ncalc\f$. The rest of the values up to
 *     \f$n=0,1,\ldots,nb-1\f$ is calculated to a lower accuracy.
 *
 * \acknowledgement
 *
 * This program is based on a program written by David J. Sookne [2] that
 * computes values of the Bessel functions \f$J_{\nu}(x)\f$ or \f$I_{\nu}(x)\f$
 * for real argument \f$x\f$ and integer order \f$\nu\f$. modifications include
 * the restriction of the computation to the Bessel function \f$I_{\nu}(x)\f$
 * for non-negative real argument, the extension of the computation to arbitrary
 * non-negative orders \f$\nu\f$, and the elimination of most underflow.
 *
 * References:
 * [1] F. W. J. Olver and D. J. Sookne, A note on backward recurrence
 *   algorithms", Math. Comput. (26), 1972, pp 125 -- 132.
 * [2] D. J. Sookne, "Bessel functions of real argument and int order", NBS
 *   Jour. of Res. B. (77B), 1973, pp. 125 -- 132.
 *
 * Modified by W. J. Cody, Applied Mathematics Division, Argonne National
 *   Laboratory, Argonne, IL, 60439, USA
 *
 * Modified by Jens Keiner, Institute of Mathematics, University of Lübeck,
 *   23560 Lübeck, Germany
 */
int nfft_smbi(const R x, const R alpha, const int nb, const int ize, R *b)
{
  /* machine dependent parameters */
  /* NSIG   - DECIMAL SIGNIFICANCE DESIRED.  SHOULD BE SET TO */
  /*          IFIX(ALOG10(2)*NBIT+1), WHERE NBIT IS THE NUMBER OF */
  /*          BITS IN THE MANTISSA OF A WORKING PRECISION VARIABLE. */
  /*          SETTING NSIG LOWER WILL RESULT IN DECREASED ACCURACY */
  /*          WHILE SETTING NSIG HIGHER WILL INCREASE CPU TIME */
  /*          WITHOUT INCREASING ACCURACY.  THE TRUNCATION ERROR */
  /*          IS LIMITED TO A RELATIVE ERROR OF T=.5*10**(-NSIG). */
  /* ENTEN  - 10.0 ** K, WHERE K IS THE LARGEST int SUCH THAT */
  /*          ENTEN IS MACHINE-REPRESENTABLE IN WORKING PRECISION. */
  /* ENSIG  - 10.0 ** NSIG. */
  /* RTNSIG - 10.0 ** (-K) FOR THE SMALLEST int K SUCH THAT */
  /*          K .GE. NSIG/4. */
  /* ENMTEN - THE SMALLEST ABS(X) SUCH THAT X/4 DOES NOT UNDERFLOW. */
  /* XLARGE - UPPER LIMIT ON THE MAGNITUDE OF X WHEN IZE=2.  BEAR */
  /*          IN MIND THAT IF ABS(X)=N, THEN AT LEAST N ITERATIONS */
  /*          OF THE BACKWARD RECURSION WILL BE EXECUTED. */
  /* EXPARG - LARGEST WORKING PRECISION ARGUMENT THAT THE LIBRARY */
  /*          EXP ROUTINE CAN HANDLE AND UPPER LIMIT ON THE */
  /*          MAGNITUDE OF X WHEN IZE=1. */
  const int nsig = MANT_DIG + 2;
  const R enten = nfft_float_property(NFFT_R_MAX);
  const R ensig = POW(K(10.0),(R)nsig);
  const R rtnsig = POW(K(10.0),-CEIL((R)nsig/K(4.0)));
  const R xlarge = K(1E4);
  const R exparg = FLOOR(LOG(POW(K(R_RADIX),K(DBL_MAX_EXP-1))));

  /* System generated locals */
  int l, n, nend, magx, nbmx, ncalc, nstart;
  R p, em, en, sum, pold, test, empal, tempa, tempb, tempc, psave, plast, tover,
    emp2al, psavel;

  magx = LRINT(FLOOR(x));

  /* return if x, nb, or ize out of range */
  if (   nb <= 0 || x < K(0.0) || alpha < K(0.0) || K(1.0) <= alpha
      || ((ize != 1 || exparg < x) && (ize != 2 || xlarge < x)))
    return (MIN(nb,0) - 1);

  /* 2-term ascending series for small x */
  if (x < rtnsig)
    return scaled_modified_bessel_i_series(x,alpha,nb,ize,b);

  ncalc = nb;
  /* forward sweep, Olver's p-sequence */

  nbmx = nb - magx;
  n = magx + 1;

  en = (R) (n+n) + (alpha+alpha);
  plast = K(1.0);
  p = en/x;

  /* significance test */
  test = ensig + ensig;

  if ((5*nsig) < (magx << 1))
    test = SQRT(test*p);
  else
    test /= POW(K(1.585),(R)magx);

  if (3 <= nbmx)
  {
    /* calculate p-sequence until n = nb-1 */
    tover = enten/ensig;
    nstart = magx+2;
    nend = nb - 1;

    for (n = nstart; n <= nend; n++)
    {
      en += K(2.0);
      pold = plast;
      plast = p;
      p = en*plast/x + pold;
      if (p > tover)
      {
        /* divide p-sequence by tover to avoid overflow. Calculate p-sequence
         * until 1 <= |p| */
        tover = enten;
        p /= tover;
        plast /= tover;
        psave = p;
        psavel = plast;
        nstart = n + 1;

        do
        {
          n++;
          en += K(2.0);
          pold = plast;
          plast = p;
          p = en*plast/x + pold;
        } while (p <= K(1.0));

        tempb = en/x;

        /* Backward test. Find ncalc as the largest n such that test is passed. */
        test = pold*plast*(K(0.5) - K(0.5)/(tempb * tempb))/ensig;
        p = plast*tover;
        n--;
        en -= K(2.0);
        nend = MIN(nb,n);

        for (ncalc = nstart; ncalc <= nend; ncalc++)
        {
          pold = psavel;
          psavel = psave;
          psave = en*psavel/x + pold;
          if (test < psave * psavel)
            break;
        }

        ncalc--;
        goto L80;
      }
    }

    n = nend;
    en = (R) (n+n) + (alpha+alpha);

    /* special significance test for 2 <= nbmx */
    test = FMAX(test,SQRT(plast*ensig)*SQRT(p+p));
  }

  /* calculate p-sequence until significance test is passed */
  do
  {
    n++;
    en += K(2.0);
    pold = plast;
    plast = p;
    p = en*plast/x + pold;
  } while (p < test);

  /* Initialize backward recursion and normalization sum. */
L80:
  n++;
  en += K(2.0);
  tempb = K(0.0);
  tempa = K(1.0)/p;
  em = (R)(n-1);
  empal = em + alpha;
  emp2al = em - K(1.0) + (alpha+alpha);
  sum = tempa*empal*emp2al/em;
  nend = n-nb;

  if (nend < 0)
  {
    /* We have n <= nb. So store b[n] and set higher orders to zero */
    b[n-1] = tempa;
    nend = -nend;
    for (l = 1; l <= nend; ++l)
      b[n-1 + l] = K(0.0);
  }
  else
  {
    if (nend != 0)
    {
      /* recur backward via difference equation, calculating b[n] until n = nb */
      for (l = 1; l <= nend; ++l)
      {
        n--;
        en -= K(2.0);
        tempc = tempb;
        tempb = tempa;
        tempa = en*tempb/x + tempc;
        em -= K(1.0);
        emp2al -= K(1.0);

        if (n == 1)
          break;

        if (n == 2)
          emp2al = K(1.0);

        empal -= K(1.0);
        sum = (sum + tempa*empal)*emp2al/em;
      }
    }

    /* store b[nb] */
    b[n-1] = tempa;

    if (nb <= 1)
    {
      sum = sum + sum + tempa;
      scaled_modified_bessel_i_normalize(x,alpha,nb,ize,b,sum);
      return ncalc;
    }

    /* calculate and store b[nb-1] */
    n--;
    en -= 2.0;
    b[n-1] = en*tempa/x + tempb;

    if (n == 1)
    {
      sum = sum + sum + b[0];
      scaled_modified_bessel_i_normalize(x,alpha,nb,ize,b,sum);
      return ncalc;
    }

    em -= K(1.0);
    emp2al -= K(1.0);

    if (n == 2)
      emp2al = K(1.0);

    empal -= K(1.0);
    sum = (sum + b[n-1]*empal)*emp2al/em;
  }

  nend = n - 2;

  if (nend != 0)
  {
    /* Calculate and store b[n] until n = 2. */
    for (l = 1; l <= nend; ++l)
    {
      n--;
      en -= K(2.0);
      b[n-1] = en*b[n]/x + b[n+1];
      em -= K(1.0);
      emp2al -= K(1.0);

      if (n == 2)
        emp2al = K(1.0);

      empal -= K(1.0);
      sum = (sum + b[n-1]*empal)*emp2al/em;
    }
  }

  /* calculate b[1] */
  b[0] = K(2.0)*empal*b[1]/x + b[2];
  sum = sum + sum + b[0];

  scaled_modified_bessel_i_normalize(x,alpha,nb,ize,b,sum);
  return ncalc;
}

int nfft_get_num_threads(void)
{
#ifdef _OPENMP
  return nfft_get_omp_num_threads();
#else
  return 1;
#endif
}

#ifdef _OPENMP
int nfft_get_omp_num_threads(void)
{
  int nthreads;
  #pragma omp parallel default(shared)
  {
    int n = omp_get_num_threads();
    #pragma omp master
    {
      nthreads = n;
    }
  }
  return nthreads;
}
#endif
