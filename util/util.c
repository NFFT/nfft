/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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
#include "nfft3.h"
#include "nfft3util.h"
#include "infft.h"

R nfft_fc(const char *cmach)
{
  const R base = FLT_RADIX;
  const R eps = EPSILON;
  const R t = MANT_DIG;
  const R emin = MIN_EXP;
  const R emax = MAX_EXP;
  const R prec = eps * base;
  static R rmin = K(1.0);
  static R rmax = K(1.0);
  const R rnd = FLTROUND;
  static R sfmin = K(-1.0);
  static short first = TRUE;

  if (first)
  {
    /* Compute rmin */
    {
      const INT n = 1 - MIN_EXP;
      INT i;
      for (i = 0; i < n; i++)
        rmin /= base;
    }
    /* Compute rmax */
    {
      INT i;
      rmax -= eps;
      for (i = 0; i < emax; i++)
        rmax *= base;
    }
    /* Compute sfmin */
    {
      R small = K(1.0) / rmax;
      sfmin = rmin;
      if (small >= sfmin)
        sfmin = small * (eps + K(1.0));
    }
    first = FALSE;
  }

  if (cmach[0] == 'E')
    return eps;
  else if (cmach[0] == 'S')
    return sfmin;
  else if (cmach[0] == 'B')
    return base;
  else if (cmach[0] == 'P')
    return prec;
  else if (cmach[0] == 'N')
    return t;
  else if (cmach[0] == 'R')
    return rnd;
  else if (cmach[0] == 'M')
    return  emin;
  else if (cmach[0] == 'U')
    return rmin;
  else if (cmach[0] == 'L')
    return emax;
  else if (cmach[0] == 'O')
    return rmax;
  else
    CK(0 /* cannot happen */);

  return K(-1.0);
} /* dlamch_ */

double nfft_elapsed_seconds(ticks t1, ticks t0)
{
  UNUSED(t1);
  UNUSED(t0);
  return elapsed(t1,t0)/TICKS_PER_SECOND;
}

int nfft_int_2_pow(int a)
{
  return (1U<< a);
}

int nfft_ld(int m)
{
  int l=0;
  int mm=m;

  while(mm>0)
    {
      mm=(mm>>1);
      l++;
    }
  return (l-1);
}

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
int nfft_next_power_of_2(int N)
{
  int n,i,logn;
  int N_is_not_power_of_2=0;

  if (N == 0)
  {
    return 1;
  }
  else
  {
    n=N;
    logn=0;
    while (n != 1)
    {
      if (n%2 == 1)
      {
          N_is_not_power_of_2=1;
      }
      n = n/2;
      logn++;
    }

    if (!N_is_not_power_of_2)
    {
      logn--;
    }

    for (i = 0; i <= logn; i++)
    {
      n = n*2;
    }

    return n;
  }
}

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
void nfft_next_power_of_2_exp(int N, int *N2, int *t)
{
  int n,i,logn;
  int N_is_not_power_of_2=0;

  if (N == 0)
  {
    *N2 = 1;
    *t = 0;
  }
  else
  {
    n=N;
    logn=0;
    while (n != 1)
    {
      if (n%2 == 1)
      {
          N_is_not_power_of_2=1;
      }
      n = n/2;
      logn++;
    }

    if (!N_is_not_power_of_2)
    {
      logn--;
    }

    for (i = 0; i <= logn; i++)
    {
      n = n*2;
    }

    *N2 = n;
    *t = logn+1;
  }
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

/** Computes integer /f$\prod_{t=0}^{d-1} v_t/f$.
 */
int nfct_prod_int(int *vec, int d)
{
  return nfft_prod_int(vec, d);
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

  sum=idx[0];
  for(t=1; t<d; t++)
    sum=sum*N[t]+idx[t];

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

double nfft_sinc(double x)
{
  if (fabs(x) < 1e-20)
    return(1.0);
  else
    return((double)(sin(x)/x));
} /* sinc */

static void bspline_help(int k, double x, double *scratch, int j, int ug, int og,
      int r)
{
  int i;                                /**< row index of the de Boor scheme  */
  int idx;                            /**< index in scratch                 */
  double a;                             /**< alpha of the de Boor scheme      */

  /* computation of one column */
  for(i=og+r-k+1, idx=og; idx>=ug; i--, idx--)
    {
      a = ((double)(x - i)) / ((double)(k - j));
      scratch[idx] = (1 - a) * scratch[idx-1] + a * scratch[idx];
    }
} /* bspline_help */

/** Computes \f$M_{k,0}\left(x\right)\f$
 *  scratch is used for de Boor's scheme
 */
double nfft_bspline(int k, double x, double *scratch)
{
  double result_value;                  /**< M_{k,0}\left(x\right)            */
  int r;                                /**< \f$x \in {\rm supp}(M_{0,r}) \f$ */
  int g1,g2;                            /**< boundaries                       */
  int j,idx,ug,og;                    /**< indices                          */
  double a;                             /**< alpha of the de Boor scheme      */

  result_value=0.0;
  if(0<x && x<k)
    {
      /* using symmetry around k/2 */
      if((k-x)<x) x=k-x;

      r=(int)(ceil(x)-1.0);

      for(idx=0; idx<k; idx++)
  scratch[idx]=0.0;

      scratch[k-r-1]=1.0;

      /* bounds of the algorithm */
      g1 = r;
      g2 = k - 1 - r;
      ug = g2;

      /* g1<=g2 holds */

      for(j=1, og=g2+1; j<=g1; j++, og++)
  {
    a = (x - r + k - 1 - og)/(k - j);
    scratch[og] = (1 - a) * scratch[og-1];
    bspline_help(k,x,scratch,j,ug+1,og-1,r);
    a = (x - r + k - 1 - ug)/(k - j);
    scratch[ug] = a * scratch[ug];
  }
      for(og-- ; j<=g2; j++)
  {
    bspline_help(k,x,scratch,j,ug+1,og,r);
    a = (x - r + k - 1 - ug)/(k - j);
    scratch[ug] = a * scratch[ug];
  }
      for( ; j<k; j++)
  {
    ug++;
    bspline_help(k,x,scratch,j,ug,og,r);
  }
      result_value = scratch[k-1];
    }
  return(result_value);
} /* bspline */

/*              mconf.h
 *
 *  Common include file for math routines
 *
 *
 *
 * SYNOPSIS:
 *
 * #include "mconf.h"
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains definitions for error codes that are
 * passed to the common error handling routine mtherr()
 * (which see).
 *
 * The file also includes a conditional assembly definition
 * for the type of computer arithmetic (IEEE, DEC, Motorola
 * IEEE, or UNKnown).
 *
 * For Digital Equipment PDP-11 and VAX computers, certain
 * IBM systems, and others that use numbers with a 56-bit
 * significand, the symbol DEC should be defined.  In this
 * mode, most floating point constants are given as arrays
 * of octal integers to eliminate decimal to binary conversion
 * errors that might be introduced by the compiler.
 *
 * For little-endian computers, such as IBM PC, that follow the
 * IEEE Standard for Binary Floating Point Arithmetic (ANSI/IEEE
 * Std 754-1985), the symbol IBMPC should be defined.  These
 * numbers have 53-bit significands.  In this mode, constants
 * are provided as arrays of hexadecimal 16 bit integers.
 *
 * Big-endian IEEE format is denoted MIEEE.  On some RISC
 * systems such as Sun SPARC, double precision constants
 * must be stored on 8-byte address boundaries.  Since integer
 * arrays may be aligned differently, the MIEEE configuration
 * may fail on such machines.
 *
 * To accommodate other types of computer arithmetic, all
 * constants are also provided in a normal decimal radix
 * which one can hope are correctly converted to a suitable
 * format by the available C language compiler.  To invoke
 * this mode, define the symbol UNK.
 *
 * An important difference among these modes is a predefined
 * set of machine arithmetic constants for each.  The numbers
 * MACHEP (the machine roundoff error), MAXNUM (largest number
 * represented), and several other parameters are preset by
 * the configuration symbol.  Check the file const.c to
 * ensure that these values are correct for your computer.
 *
 * Configurations NANS, INFINITIES, MINUSZERO, and DENORMAL
 * may fail on many systems.  Verify that they are supposed
 * to work on your computer.
 */
/*
Cephes Math Library Release 2.3:  June, 1995
Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
*/

/* Define if the `long double' type works.  */
#define HAVE_LONG_DOUBLE 1

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* Define if floating point words are bigendian.  */
/* #undef FLOAT_WORDS_BIGENDIAN */

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1

/* Name of package */
//#define PACKAGE "cephes"

/* Version number of package */
//#define VERSION "2.7"

/* Constant definitions for math error conditions
 */

#define DOMAIN    1  /* argument domain error */
#define SING    2  /* argument singularity */
#define OVERFLOW  3  /* overflow range error */
#define UNDERFLOW  4  /* underflow range error */
#define TLOSS    5  /* total loss of precision */
#define PLOSS    6  /* partial loss of precision */

#define EDOM    33
#define ERANGE    34

/* Type of computer arithmetic */

/* PDP-11, Pro350, VAX:
 */
/* #define DEC 1 */

/* Intel IEEE, low order words come first:
 */
/* #define IBMPC 1 */

/* Motorola IEEE, high order words come first
 * (Sun 680x0 workstation):
 */
/* #define MIEEE 1 */

/* UNKnown arithmetic, invokes coefficients given in
 * normal decimal format.  Beware of range boundary
 * problems (MACHEP, MAXLOG, etc. in const.c) and
 * roundoff problems in pow.c:
 * (Sun SPARCstation)
 */
#define UNK 1

/* If you define UNK, then be sure to set BIGENDIAN properly. */
#ifdef FLOAT_WORDS_BIGENDIAN
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif
/* Define this `volatile' if your compiler thinks
 * that floating point arithmetic obeys the associative
 * and distributive laws.  It will defeat some optimizations
 * (but probably not enough of them).
 *
 * #define VOLATILE volatile
 */
#define VOLATILE

/* For 12-byte long doubles on an i386, pad a 16-bit short 0
 * to the end of real constants initialized by integer arrays.
 *
 * #define XPD 0,
 *
 * Otherwise, the type is 10 bytes long and XPD should be
 * defined blank (e.g., Microsoft C).
 *
 * #define XPD
 */
#define XPD 0,

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to ask for infinity support, else undefine. */
#define INFINITIES 1

/* Define to ask for support of numbers that are Not-a-Number,
   else undefine.  This may automatically define INFINITIES in some files. */
#define NANS 1

/* Define to distinguish between -0.0 and +0.0.  */
#define MINUSZERO 1

/* Define 1 for ANSI C atan2() function
   See atan.c and clog.c. */
#define ANSIC 1

/*              chbevl.c
 *
 *  Evaluate Chebyshev series
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N], chebevl();
 *
 * y = chbevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the series
 *
 *        N-1
 *         - '
 *  y  =   >   coef[i] T (x/2)
 *         -            i
 *        i=0
 *
 * of Chebyshev polynomials Ti at argument x/2.
 *
 * Coefficients are stored in reverse order, i.e. the zero
 * order term is last in the array.  Note N is the number of
 * coefficients, not the order.
 *
 * If coefficients are for the interval a to b, x must
 * have been transformed to x -> 2(2x - b - a)/(b-a) before
 * entering the routine.  This maps x from (a, b) to (-1, 1),
 * over which the Chebyshev polynomials are defined.
 *
 * If the coefficients are for the inverted interval, in
 * which (a, b) is mapped to (1/b, 1/a), the transformation
 * required is x -> 2(2ab/x - b - a)/(b-a).  If b is infinity,
 * this becomes x -> 4a/x - 1.
 *
 *
 *
 * SPEED:
 *
 * Taking advantage of the recurrence properties of the
 * Chebyshev polynomials, the routine requires one more
 * addition per loop than evaluating a nested polynomial of
 * the same degree.
 *
 */

/*              chbevl.c  */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1985, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

static double chbevl(double x, double *array, int n)
{
  double b0, b1, b2, *p;
  int i;

  p = array;
  b0 = *p++;
  b1 = 0.0;
  i = n - 1;

  do
    {
      b2 = b1;
      b1 = b0;
      b0 = x * b1  -  b2  + *p++;
    }
  while( --i );

  return( 0.5*(b0-b2) );
}

/*              i0.c
 *
 *  Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0();
 *
 * y = i0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         6000       8.2e-17     1.9e-17
 *    IEEE      0,30        30000       5.8e-16     1.4e-16
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */

#ifdef UNK
static double A[] =
{
-4.41534164647933937950E-18,
 3.33079451882223809783E-17,
-2.43127984654795469359E-16,
 1.71539128555513303061E-15,
-1.16853328779934516808E-14,
 7.67618549860493561688E-14,
-4.85644678311192946090E-13,
 2.95505266312963983461E-12,
-1.72682629144155570723E-11,
 9.67580903537323691224E-11,
-5.18979560163526290666E-10,
 2.65982372468238665035E-9,
-1.30002500998624804212E-8,
 6.04699502254191894932E-8,
-2.67079385394061173391E-7,
 1.11738753912010371815E-6,
-4.41673835845875056359E-6,
 1.64484480707288970893E-5,
-5.75419501008210370398E-5,
 1.88502885095841655729E-4,
-5.76375574538582365885E-4,
 1.63947561694133579842E-3,
-4.32430999505057594430E-3,
 1.05464603945949983183E-2,
-2.37374148058994688156E-2,
 4.93052842396707084878E-2,
-9.49010970480476444210E-2,
 1.71620901522208775349E-1,
-3.04682672343198398683E-1,
 6.76795274409476084995E-1
};
#endif

#ifdef DEC
static unsigned short A[] = {
0121642,0162671,0004646,0103567,
0022431,0115424,0135755,0026104,
0123214,0023533,0110365,0156635,
0023767,0033304,0117662,0172716,
0124522,0100426,0012277,0157531,
0025254,0155062,0054461,0030465,
0126010,0131143,0013560,0153604,
0026517,0170577,0006336,0114437,
0127227,0162253,0152243,0052734,
0027724,0142766,0061641,0160200,
0130416,0123760,0116564,0125262,
0031066,0144035,0021246,0054641,
0131537,0053664,0060131,0102530,
0032201,0155664,0165153,0020652,
0132617,0061434,0074423,0176145,
0033225,0174444,0136147,0122542,
0133624,0031576,0056453,0020470,
0034211,0175305,0172321,0041314,
0134561,0054462,0147040,0165315,
0035105,0124333,0120203,0162532,
0135427,0013750,0174257,0055221,
0035726,0161654,0050220,0100162,
0136215,0131361,0000325,0041110,
0036454,0145417,0117357,0017352,
0136702,0072367,0104415,0133574,
0037111,0172126,0072505,0014544,
0137302,0055601,0120550,0033523,
0037457,0136543,0136544,0043002,
0137633,0177536,0001276,0066150,
0040055,0041164,0100655,0010521
};
#endif

#ifdef IBMPC
static unsigned short A[] = {
0xd0ef,0x2134,0x5cb7,0xbc54,
0xa589,0x977d,0x3362,0x3c83,
0xbbb4,0x721e,0x84eb,0xbcb1,
0x5eba,0x93f6,0xe6d8,0x3cde,
0xfbeb,0xc297,0x5022,0xbd0a,
0x2627,0x4b26,0x9b46,0x3d35,
0x1af0,0x62ee,0x164c,0xbd61,
0xd324,0xe19b,0xfe2f,0x3d89,
0x6abc,0x7a94,0xfc95,0xbdb2,
0x3c10,0xcc74,0x98be,0x3dda,
0x9556,0x13ae,0xd4fe,0xbe01,
0xcb34,0xa454,0xd903,0x3e26,
0x30ab,0x8c0b,0xeaf6,0xbe4b,
0x6435,0x9d4d,0x3b76,0x3e70,
0x7f8d,0x8f22,0xec63,0xbe91,
0xf4ac,0x978c,0xbf24,0x3eb2,
0x6427,0xcba5,0x866f,0xbed2,
0x2859,0xbe9a,0x3f58,0x3ef1,
0x1d5a,0x59c4,0x2b26,0xbf0e,
0x7cab,0x7410,0xb51b,0x3f28,
0xeb52,0x1f15,0xe2fd,0xbf42,
0x100e,0x8a12,0xdc75,0x3f5a,
0xa849,0x201a,0xb65e,0xbf71,
0xe3dd,0xf3dd,0x9961,0x3f85,
0xb6f0,0xf121,0x4e9e,0xbf98,
0xa32d,0xcea8,0x3e8a,0x3fa9,
0x06ea,0x342d,0x4b70,0xbfb8,
0x88c0,0x77ac,0xf7ac,0x3fc5,
0xcd8d,0xc057,0x7feb,0xbfd3,
0xa22a,0x9035,0xa84e,0x3fe5,
};
#endif

#ifdef MIEEE
static unsigned short A[] = {
0xbc54,0x5cb7,0x2134,0xd0ef,
0x3c83,0x3362,0x977d,0xa589,
0xbcb1,0x84eb,0x721e,0xbbb4,
0x3cde,0xe6d8,0x93f6,0x5eba,
0xbd0a,0x5022,0xc297,0xfbeb,
0x3d35,0x9b46,0x4b26,0x2627,
0xbd61,0x164c,0x62ee,0x1af0,
0x3d89,0xfe2f,0xe19b,0xd324,
0xbdb2,0xfc95,0x7a94,0x6abc,
0x3dda,0x98be,0xcc74,0x3c10,
0xbe01,0xd4fe,0x13ae,0x9556,
0x3e26,0xd903,0xa454,0xcb34,
0xbe4b,0xeaf6,0x8c0b,0x30ab,
0x3e70,0x3b76,0x9d4d,0x6435,
0xbe91,0xec63,0x8f22,0x7f8d,
0x3eb2,0xbf24,0x978c,0xf4ac,
0xbed2,0x866f,0xcba5,0x6427,
0x3ef1,0x3f58,0xbe9a,0x2859,
0xbf0e,0x2b26,0x59c4,0x1d5a,
0x3f28,0xb51b,0x7410,0x7cab,
0xbf42,0xe2fd,0x1f15,0xeb52,
0x3f5a,0xdc75,0x8a12,0x100e,
0xbf71,0xb65e,0x201a,0xa849,
0x3f85,0x9961,0xf3dd,0xe3dd,
0xbf98,0x4e9e,0xf121,0xb6f0,
0x3fa9,0x3e8a,0xcea8,0xa32d,
0xbfb8,0x4b70,0x342d,0x06ea,
0x3fc5,0xf7ac,0x77ac,0x88c0,
0xbfd3,0x7feb,0xc057,0xcd8d,
0x3fe5,0xa84e,0x9035,0xa22a
};
#endif


/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

#ifdef UNK
static double B[] =
{
-7.23318048787475395456E-18,
-4.83050448594418207126E-18,
 4.46562142029675999901E-17,
 3.46122286769746109310E-17,
-2.82762398051658348494E-16,
-3.42548561967721913462E-16,
 1.77256013305652638360E-15,
 3.81168066935262242075E-15,
-9.55484669882830764870E-15,
-4.15056934728722208663E-14,
 1.54008621752140982691E-14,
 3.85277838274214270114E-13,
 7.18012445138366623367E-13,
-1.79417853150680611778E-12,
-1.32158118404477131188E-11,
-3.14991652796324136454E-11,
 1.18891471078464383424E-11,
 4.94060238822496958910E-10,
 3.39623202570838634515E-9,
 2.26666899049817806459E-8,
 2.04891858946906374183E-7,
 2.89137052083475648297E-6,
 6.88975834691682398426E-5,
 3.36911647825569408990E-3,
 8.04490411014108831608E-1
};
#endif

#ifdef DEC
static unsigned short B[] = {
0122005,0066672,0123124,0054311,
0121662,0033323,0030214,0104602,
0022515,0170300,0113314,0020413,
0022437,0117350,0035402,0007146,
0123243,0000135,0057220,0177435,
0123305,0073476,0144106,0170702,
0023777,0071755,0017527,0154373,
0024211,0052214,0102247,0033270,
0124454,0017763,0171453,0012322,
0125072,0166316,0075505,0154616,
0024612,0133770,0065376,0025045,
0025730,0162143,0056036,0001632,
0026112,0015077,0150464,0063542,
0126374,0101030,0014274,0065457,
0127150,0077271,0125763,0157617,
0127412,0104350,0040713,0120445,
0027121,0023765,0057500,0001165,
0030407,0147146,0003643,0075644,
0031151,0061445,0044422,0156065,
0031702,0132224,0003266,0125551,
0032534,0000076,0147153,0005555,
0033502,0004536,0004016,0026055,
0034620,0076433,0142314,0171215,
0036134,0146145,0013454,0101104,
0040115,0171425,0062500,0047133
};
#endif

#ifdef IBMPC
static unsigned short B[] = {
0x8b19,0x54ca,0xadb7,0xbc60,
0x9130,0x6611,0x46da,0xbc56,
0x8421,0x12d9,0xbe18,0x3c89,
0x41cd,0x0760,0xf3dd,0x3c83,
0x1fe4,0xabd2,0x600b,0xbcb4,
0xde38,0xd908,0xaee7,0xbcb8,
0xfb1f,0xa3ea,0xee7d,0x3cdf,
0xe6d7,0x9094,0x2a91,0x3cf1,
0x629a,0x7e65,0x83fe,0xbd05,
0xbb32,0xcf68,0x5d99,0xbd27,
0xc545,0x0d5f,0x56ff,0x3d11,
0xc073,0x6b83,0x1c8c,0x3d5b,
0x8cec,0xfa26,0x4347,0x3d69,
0x8d66,0x0317,0x9043,0xbd7f,
0x7bf2,0x357e,0x0fd7,0xbdad,
0x7425,0x0839,0x511d,0xbdc1,
0x004f,0xabe8,0x24fe,0x3daa,
0x6f75,0xc0f4,0xf9cc,0x3e00,
0x5b87,0xa922,0x2c64,0x3e2d,
0xd56d,0x80d6,0x5692,0x3e58,
0x616e,0xd9cd,0x8007,0x3e8b,
0xc586,0xc101,0x412b,0x3ec8,
0x9e52,0x7899,0x0fa3,0x3f12,
0x9049,0xa2e5,0x998c,0x3f6b,
0x09cb,0xaca8,0xbe62,0x3fe9
};
#endif

#ifdef MIEEE
static unsigned short B[] = {
0xbc60,0xadb7,0x54ca,0x8b19,
0xbc56,0x46da,0x6611,0x9130,
0x3c89,0xbe18,0x12d9,0x8421,
0x3c83,0xf3dd,0x0760,0x41cd,
0xbcb4,0x600b,0xabd2,0x1fe4,
0xbcb8,0xaee7,0xd908,0xde38,
0x3cdf,0xee7d,0xa3ea,0xfb1f,
0x3cf1,0x2a91,0x9094,0xe6d7,
0xbd05,0x83fe,0x7e65,0x629a,
0xbd27,0x5d99,0xcf68,0xbb32,
0x3d11,0x56ff,0x0d5f,0xc545,
0x3d5b,0x1c8c,0x6b83,0xc073,
0x3d69,0x4347,0xfa26,0x8cec,
0xbd7f,0x9043,0x0317,0x8d66,
0xbdad,0x0fd7,0x357e,0x7bf2,
0xbdc1,0x511d,0x0839,0x7425,
0x3daa,0x24fe,0xabe8,0x004f,
0x3e00,0xf9cc,0xc0f4,0x6f75,
0x3e2d,0x2c64,0xa922,0x5b87,
0x3e58,0x5692,0x80d6,0xd56d,
0x3e8b,0x8007,0xd9cd,0x616e,
0x3ec8,0x412b,0xc101,0xc586,
0x3f12,0x0fa3,0x7899,0x9e52,
0x3f6b,0x998c,0xa2e5,0x9049,
0x3fe9,0xbe62,0xaca8,0x09cb
};
#endif

/** Modified Bessel function of order zero.
 *  Cephes Math Library Release 2.8:  June, 2000
 *  Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */
double nfft_i0(double x)
{
  double y;

  if( x < 0 )
    x = -x;
  if( x <= 8.0 )
    {
      y = (x/2.0) - 2.0;
      return( exp(x) * chbevl( y, A, 30 ) );
    }

  return(  exp(x) * chbevl( 32.0/x - 2.0, B, 25 ) / sqrt(x) );
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

static double l_1_complex(double _Complex *x, double _Complex *y, int n)
{
  int k;
  double l1;

  if(y==NULL)
    for(k=0,l1=0; k<n; k++)
      l1+=cabs(x[k]);
  else
    for(k=0,l1=0; k<n; k++)
      l1+=cabs(x[k]-y[k]);

  return l1;
}

static double l_1_double(double *x, double *y, int n)
{
  int k;
  double l1;

  if(y==NULL)
    for(k=0,l1=0; k<n; k++)
      l1+=fabs(x[k]);
  else
    for(k=0,l1=0; k<n; k++)
      l1+=fabs(x[k]-y[k]);

  return l1;
}




static double l_2_complex(double _Complex *x, double _Complex *y, int n)
{
  int k;
  double l22;

  if(y==NULL)
    for(k=0,l22=0; k<n; k++)
      l22+=conj(x[k])*x[k];
  else
    for(k=0,l22=0; k<n; k++)
      l22+=conj(x[k]-y[k])*(x[k]-y[k]);

  return sqrt(l22);
}

static double l_2_double(double *x, double *y, int n)
{
  int k;
  double l22;

  if(y==NULL)
    for(k=0,l22=0; k<n; k++)
      l22+=x[k]*x[k];
  else
    for(k=0,l22=0; k<n; k++)
      l22+=(x[k]-y[k])*(x[k]-y[k]);

  return sqrt(l22);
}

static double l_infty_complex(double _Complex *x, double _Complex *y, int n)
{
  int k;
  double linfty;

  if(y==NULL)
    for(k=0,linfty=0; k<n; k++)
      linfty=((linfty<cabs(x[k]))?cabs(x[k]):linfty);
  else
    for(k=0,linfty=0; k<n; k++)
      linfty=((linfty<cabs(x[k]-y[k]))?cabs(x[k]-y[k]):linfty);

  return linfty;
}


static double l_infty_double(double *x, double *y, int n)
{
  int k;
  double linfty;

  if(y==NULL)
    for(k=0,linfty=0; k<n; k++)
      linfty=((linfty<fabs(x[k]))?fabs(x[k]):linfty);
  else
    for(k=0,linfty=0; k<n; k++)
      linfty=((linfty<fabs(x[k]-y[k]))?fabs(x[k]-y[k]):linfty);

  return linfty;
}





/** computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
double nfft_error_l_infty_complex(double _Complex *x, double _Complex *y, int n)
{
  return (l_infty_complex(x, y, n)/l_infty_complex(x, NULL, n));
}

/** computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
double nfft_error_l_infty_double(double *x, double *y, int n)
{
  return (l_infty_double(x, y, n)/l_infty_double(x, NULL, n));
}



/** computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
double nfft_error_l_infty_1_complex(double _Complex *x, double _Complex *y, int n,
             double _Complex *z, int m)
{
  return (l_infty_complex(x, y, n)/l_1_complex(z, NULL, m));
}

/** computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
double nfft_error_l_infty_1_double(double *x, double *y, int n,
                              double *z, int m)
{
  return (l_infty_double(x, y, n)/l_1_double(z, NULL, m));
}



/** computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
double nfft_error_l_2_complex(double _Complex *x, double _Complex *y, int n)
{
  return (l_2_complex(x, y, n)/l_2_complex(x, NULL, n));
}

/** computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
double  nfft_error_l_2_double(double *x, double *y, int n)
{
  return (l_2_double(x, y, n)/l_2_double(x, NULL, n));
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

void X(vrand_unit_complex)(C *x, const int n)
{
  int k;

  for (k = 0; k < n; k++)
    x[k] = nfft_drand48() + II*nfft_drand48();
}

void X(vrand_shifted_unit_double)(R *x, const int n)
{
  int k;

  for (k = 0; k < n; k++)
    x[k] = nfft_drand48() - K(0.5);
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
  const R enmten = K(4.0)*nfft_fc("U");
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
  const R enmten = K(4.0)*nfft_fc("U");
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
  const R enten = nfft_fc("O");//POW(K(10.0),K(R_MAX_10_EXP));
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

/* Coefficients for Lanzcos's approximation to the Gamma function. Can be
 * regenerated with Mathematica from file lambda.nb. */
#if defined(NFFT_LDOUBLE)
  #if SIZEOF_LONG_DOUBLE == 16
    static const R lanzcos[22] =
    {
        K(6.992159592766260467081279017133763989752796675052e-10),
        K(2.8123669719981441131615873659455441181381798148252),
        K(-20.710427744802709693741006921011198460802137523364),
        K(68.91604215552677659752923125507927632926035487156),
        K(-137.06264623052621370808270752038118534584695692919),
        K(181.60793022155827236320731899262871987253012010211),
        K(-169.18458743339031548467251800473390300200144720615),
        K(114.00119644872831703640105977147389190644587171368),
        K(-56.31458215179027632548005714718150178553786189071),
        K(20.455656046361896023848762777393331001159409886467),
        K(-5.4335605693344340005752804679314670882396415124971),
        K(1.0410515397995220235480840944608088445246607731035),
        K(-0.14064836221797762435665594700765635768382988137516),
        K(0.01295767486030631796076454893320245106478793545067),
        K(-0.00077606781904462836224565012148234163222031403439914),
        K(0.000028229941587194672222791862511243708412223569786751),
        K(-5.6505789508859496251539126911727595443972191305184e-7),
        K(5.3642107416477670941393734212506122918693783684933e-9),
        K(-1.9054477424755607958337633536133218546180235027896e-11),
        K(1.6767560615458281885574906761496245907479292929211e-14),
        K(-1.6069159098936011461868918827916292610798729856278e-18),
        K(2.092833837734494231412616008353349620313749127668e-24)};
    static const R g = K(22.);
  #elif SIZEOF_LONG_DOUBLE == 12
  static const R lanzcos[16] =
  {
    K(0.9999999999999999999998151990326234090755978504801206971875240538249),
    K(42298.205314163579685884559965825876032898966930749727725308255),
    K(-139068.033365719940583314330900395851082150878479801707616019365),
    K(181526.445910064425982455538776735270618298440139068097226843822),
    K(-119998.518895134887942904852888927609503242910142130425166120992),
    K(42498.279665085565336840634629253119065594761580523677606897758),
    K(-7861.992932117389960417642813146578778617649628454257358410924),
    K(687.805169061068980481146382924730399157260064751317893043547),
    K(-23.0082756959823209227391335933370375125546087097458511646070585787792),
    K(0.1864230791562731665361253542681001772953257132731654277546996312946),
    K(-0.0001094217900494200137669846537980382942701406476771516598942134305),
    K(3.8119733228437385220001515981932509556399175176128943412654167182529e-10),
    K(-4.7545060695949573669792060004330071744683150710413248650156210893946e-10),
    K(2.3162039672883094408744543321924181587375284941017905602597030297461e-10),
    K(-5.469699080479773550184561017808390160175857819276162113300096555236e-11),
    K(5.761353296381615770264101368671565767532131044262156102602431704e-12)
  };
  static const R g = K(10.900511);
  #else
    #error Unsupported size of long double
  #endif
#elif defined(NFFT_SINGLE)
  static const R lanzcos[5] =
  {
    0.032648921146387503785008589380436838603631779746936L,
    1.1886889953894899130227639793483334640624981688106L,
    -1.0684103609197946791553357282032624745327225621591L,
    0.18871201360609388869017219868989358810366207478127L,
    -0.002745022750941733257487411560167743487986763037079L
  };
  static const R g = K(4.3408820000000000000000000000000000000000000000000);
#else
  static const R lanzcos[13] = {
  K(0.006061842346248906525783753964555936883222463665497),
  K(1.4256283548148038891120508241039168186977000466332),
  K(-2.1475341954959013287000186097849561222515097674559),
  K(0.95726691187623991672746418765475493079905679642641),
  K(-0.12868865928482566615781458898635486929956555429179),
  K(0.0030886434783642720644480000315021785318956671151758),
  K(-9.7849610972063267859543785146578805638938457810257e-7),
  K(-1.7768848405624834159045007110755286743551855068815e-8),
  K(1.9851561978478612433621720534889693635174319844665e-8),
  K(-1.2477453878744212747015501336682760715977546596628e-8),
  K(5.6096232232009359480838882000975403703116221325714e-9),
  K(-1.6133242708801833407734306938668413119166269573061e-9),
  K(2.191097792670207503670965523347118763341507005123e-10)};
  static const R lanzcosp[12] = {
  K(3.0627152859269464002981250710810845422307962730128e8),
  K(4.775603155294306744694136459666282665368397175884e8),
  K(3.3795733209388037978795882754036968577436873787256e8),
  K(1.4313371466803237581908862806209047500008516090763e8),
  K(4.026600756736916571780374444713775549949468190092e7),
  K(7.8910786118305655779880429845306886138919766854132e6),
  K(1.0980346455035079893427073585191634273991944222518e6),
  K(108372.57068540043287489056427702756001287882054205),
  K(7427.8953058166394494408791055208785413418341223703),
  K(336.44978500306352175687290720400758052504709702257),
  K(9.0582421437926674355635744116038966456459610306128),
  K(0.10976007071323978811079010281977761670657790941656)};
  static const R g = K(6.0246800407767295837402343750000000000000000000000);
#endif

/* Computes the partial fractions sum in Lanzcos' approximation. */
static inline R csum(const R z)
{
  const short n = sizeof(lanzcos)/sizeof(lanzcos[0]);
  R r = K(0.0), y = z + (R)(n-1);
  short i;
  for (i = n-1; i >= 1; i--, y -= K(1.0))
    r += lanzcos[i]/y;
  r += lanzcos[0];
  return r;
}

static inline R csump(const R z)
{
  const short n = sizeof(lanzcos)/sizeof(lanzcos[0]);
  R num = K(0.0), denom = K(1.0);
  short i;
  for (i = n-1; i >= 1; i--)
  {
    num = z * num + lanzcosp[i-1];
    denom *= z + (R)(i);
  }
  return lanzcos[0] + (num/denom);
}

/* Euler's number */
static DK(KE,2.7182818284590452353602874713526624977572470937000);

R nfft_lambda(const R z, const R eps)
{
  R d = K(1.0) - eps, zpg = z + g, emh = eps - K(0.5);
  return EXP(-LOG1P(d/(zpg+emh))*(z+emh)) * POW(KE/(zpg+K(0.5)),d)
    * (csump(z-d)/csump(z));
}

/* Computes lambda2(mu, nu) = Sqrt(Gamma(mu+nu+1)/(Gamma(mu+1)Gamma(nu+1)))
 * using Lanczos' approximation. */
R nfft_lambda2(const R mu, const R nu)
{
  if (mu == K(0.0) || nu == K(0.0))
    return K(1.0);
  else
    return
      SQRT(
        POW((mu+nu+g+K(0.5))/(K(1.0)*(mu+g+K(0.5))),mu)
      * POW((mu+nu+g+K(0.5))/(K(1.0)*(nu+g+K(0.5))),nu)
      * SQRT(KE*(mu+nu+g+K(0.5))/((mu+g+K(0.5))*(nu+g+K(0.5))))
      * (csum(mu+nu)/(csum(mu)*csum(nu)))
      );
}

/**
 * Prints an error message for a failed assertion together with filename and the
 * line where the assertion failed.
 */
void nfft_assertion_failed(const char *s, int line, const char *file)
{
  fflush(stdout);
  fprintf(stderr, "fmmw: %s:%d: assertion failed: %s\n", file, line, s);
#ifdef HAVE_ABORT
  /* Use abort function. */
  abort();
#else
  /* Use exit function. */
  exit(EXIT_FAILURE);
#endif
}

/* We declare drand48() and srand48() ourselves, if they are is not declared in
 * math.h (e.g. on SuSE 9.3), grrr. */
#include "config.h"
#if HAVE_DECL_DRAND48 == 0
  extern double drand48(void);
#endif
#if HAVE_DECL_SRAND48 == 0
  extern void srand48(long int);
#endif

double nfft_drand48(void)
{
#ifdef HAVE_DRAND48
  return drand48();
#else
  return ((R)rand())/((R)RAND_MAX);
#endif
}

void nfft_srand48(long int seed)
{
#ifdef HAVE_SRAND48
  srand48(seed);
#else
  srand((unsigned int)seed);
#endif
}
