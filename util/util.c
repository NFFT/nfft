/** Sources for utilities.
 *  functions for vectors, window functions, ...
 *  (c) if not stated otherwise: Daniel Potts, Stefan Kunis
 */

#include "config.h"

#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef HAVE_MALLINFO
#  include <malloc.h>
#endif
//#include <time.h>


/** Actual used CPU time in seconds.
 *  Calls getrusage, limited accuracy
 */
double second()
{
  static struct rusage temp;
  double foo, foo1;
  
  getrusage(RUSAGE_SELF,&temp);
  foo     = temp.ru_utime.tv_sec;       /* seconds                           */
  foo1    = temp.ru_utime.tv_usec;      /* uSecs                             */
  return  foo  + (foo1/1000000.0);      /* milliseconds                      */
}

#ifdef HAVE_TOTAL_USED_MEMORY
int total_used_memory()
{
  struct mallinfo m;
  m=mallinfo();
  return m.hblkhd + m.uordblks;
}
#endif

int int_2_pow(int a)
{
  return (1U<< a);
}

int ld(int m)
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
int next_power_of_2(int N)
{
  int n,i,logn;
  int N_is_not_power_of_2=0;

  n=N;
  logn=0;
  while(n!=1)
    {
      if(n%2==1)
	N_is_not_power_of_2=1;
      n=n/2;
      logn++;
    }

  if(!N_is_not_power_of_2)
    logn--;

  for(i=0;i<=logn;i++)
    n=n*2;

  return n;
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

/** Sinus cardinalis
 *  \f$\frac{sin\left(x\right)}{x}$ 
 */
double sinc(double x)
{
  if (fabs(x) < 1e-20)
    return(1.0);
  else
    return((double)(sin(x)/x));
} /* sinc */

void bspline_help(int k, double x, double *scratch, int j, int ug, int og,
		  int r)
{
  int i;                                /**< row index of the de Boor scheme  */
  int index;                            /**< index in scratch                 */
  double a;                             /**< alpha of the de Boor scheme      */
  
  /* computation of one column */
  for(i=og+r-k+1, index=og; index>=ug; i--, index--)
    {
      a = ((double)(x - i)) / ((double)(k - j));
      scratch[index] = (1 - a) * scratch[index-1] + a * scratch[index];
    }
} /* bspline_help */

/** Computes \f$M_{k,0}\left(x\right)\f$
 *  scratch is used for de Boor's scheme
 */
double bspline(int k, double x, double *scratch)
{
  double result_value;                  /**< M_{k,0}\left(x\right)            */
  int r;                                /**< \f$x \in {\rm supp}(M_{0,r}) \f$ */
  int g1,g2;                            /**< boundaries                       */ 
  int j,index,ug,og;                    /**< indices                          */
  double a;                             /**< alpha of the de Boor scheme      */

  result_value=0.0;
  if(0<x && x<k)
    {
      /* using symmetry around k/2 */
      if((k-x)<x) x=k-x;

      r=(int)(ceil(x)-1.0);
      
      for(index=0; index<k; index++)
	scratch[index]=0.0;
	
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

/*							mconf.h
 *
 *	Common include file for math routines
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
#define PACKAGE "cephes"

/* Version number of package */
#define VERSION "2.7"

/* Constant definitions for math error conditions
 */

#define DOMAIN		1	/* argument domain error */
#define SING		2	/* argument singularity */
#define OVERFLOW	3	/* overflow range error */
#define UNDERFLOW	4	/* underflow range error */
#define TLOSS		5	/* total loss of precision */
#define PLOSS		6	/* partial loss of precision */

#define EDOM		33
#define ERANGE		34

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

/*							chbevl.c
 *
 *	Evaluate Chebyshev series
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

/*							chbevl.c	*/

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1985, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

double chbevl(double x, double array[], int n)
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

/*							i0.c
 *
 *	Modified Bessel function of order zero
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
double i0(double x)
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
double dot_complex(complex *x, int n)
{
  int k;
  double dot;

  for(k=0,dot=0; k<n; k++)
    dot+=conj(x[k])*x[k];

  return dot;
}

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double dot_w_complex(complex *x, double *w, int n)
{
  int k;
  double dot;

  for(k=0,dot=0.0; k<n; k++)
    dot+=w[k]*conj(x[k])*x[k];

  return dot;
}

/** Computes the weighted inner/dot product 
    \f$x^H (w\odot w2\odot w2 \odot x)\f$.
 */
double dot_w_w2_complex(complex *x, double *w, double *w2, int n)
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
double dot_w2_complex(complex *x, double *w2, int n)
{
  int k;
  double dot;

  for(k=0,dot=0.0; k<n; k++)
    dot+=w2[k]*w2[k]*conj(x[k])*x[k];

  return dot;
}

/** Copies \f$x \leftarrow y\f$.
 */
void cp_complex(complex *x, complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=y[k];
}

/** Copies \f$x \leftarrow y\f$.
 */
void cp_double(double *x, double *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=y[k];
}

/** Copies \f$x \leftarrow a y\f$.
 */
void cp_a_complex(complex *x, double a, complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*y[k];
}

/** Copies \f$x \leftarrow w\odot y\f$.
 */
void cp_w_complex(complex *x, double *w, complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=w[k]*y[k];
}

/** Updates \f$x \leftarrow a x + y\f$.
 */
void upd_axpy_complex(complex *x, double a, complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+y[k];
}

/** Updates \f$x \leftarrow x + a y\f$.
 */
void upd_xpay_complex(complex *x, double a, complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]+=a*y[k];
}

/** Updates \f$x \leftarrow a x + b y\f$.
 */
void upd_axpby_complex(complex *x, double a, complex *y, double b, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+b*y[k];
}

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void upd_xpawy_complex(complex *x, double a, double *w, complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]+=a*w[k]*y[k];
}

/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void upd_axpwy_complex(complex *x, double a, double *w, complex *y, int n)
{
  int k;

  for(k=0;k<n;k++)
    x[k]=a*x[k]+w[k]*y[k];
}

void fftshift_complex(complex *x, int d, int* N)
{
  int d_pre, d_act, d_post;
  int N_pre, N_act, N_post;
  int k_pre, k_act, k_post;
  int k,k_swap;

  complex x_swap;

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

double l_1_complex(complex *x, complex *y, int n)
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

double l_2_complex(complex *x, complex *y, int n)
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

double l_infty_complex(complex *x, complex *y, int n)
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

/** computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
double error_l_infty_complex(complex *x, complex *y, int n)
{
  return (l_infty_complex(x, y, n)/l_infty_complex(x, NULL, n));
}

/** computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
double error_l_infty_1_complex(complex *x, complex *y, int n,
			       complex *z, int m)
{
  return (l_infty_complex(x, y, n)/l_1_complex(z, NULL, m));
}

/** computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
double error_l_2_complex(complex *x, complex *y, int n)
{
  return (l_2_complex(x, y, n)/l_2_complex(x, NULL, n));
}

/** vector print
 */
void vpr_int(int *x, int n, char *text)
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

/** vector print
 */
void vpr_double(double *x, int n, char *text)
{
  int k;

  if(x==NULL)
    {
      printf("null pointer\n");
      fflush(stdout);
      exit(-1);
    }

  if(text!=NULL)
  {
      printf ("\n %s, adr=%p\n", text, (void*)x);
      for (k=0; k<n; k++)
      {
	  if (k%8==0) 
	      printf("%6d.\t", k);
	  printf("%+.1E,", x[k]);
	  if (k%8==7) 
	      printf("\n");
      }
      if (n%8!=0)
        printf("\n");
  }
  else
      for (k=0; k<n; k++)
	  printf("%+E,\n", x[k]);
  fflush(stdout);
}

/** vector print
 */
void vpr_complex(complex *x, int n, char *text)
{
  int k;

  if(text!=NULL)
  {
      printf ("\n %s, adr=%p\n", text, (void*)x);
      for (k=0; k<n; k++)
      {
	  if (k%4==0) 
	      printf("%6d.\t", k);
	  printf("%+.1E%+.1Ei,", creal(x[k]), cimag(x[k]));
	  if (k%4==3) 
	      printf("\n");
      }
      if (n%4!=0)
        printf("\n");
  }
  else
      for (k=0; k<n; k++)
	  printf("%+E%+Ei,\n", creal(x[k]), cimag(x[k]));
  fflush(stdout);
}

void vrand_unit_complex(complex *x, int n)
{
  int k;

  for (k=0; k<n; k++)
    x[k] = ((double)rand())/RAND_MAX + I*((double)rand())/RAND_MAX;
}

void vrand_shifted_unit_double(double *x, int n)
{
  int k;

  for (k=0; k<n; k++)
    x[k] = ((double)rand())/RAND_MAX - 0.5;
}


/** Computes non periodic voronoi weights 
 *  assumes ordered x_j */
void voronoi_weights_1d(double *w, double *x, int M)
{
  int j;
  
  w[0]=(x[1]-x[0])/2;
  for(j=1;j<M-1;j++)
    w[j]=(x[j+1]-x[j-1])/2;
  w[M-1]=(x[M-1]-x[M-2])/2;
}

/** Computes the damping factor for the modified Fejer kernel.
 *  /f$\frac{2}{N}\left(1-\frac{\left|2k+1\right|}{N}\right)/f$
 */
double modified_fejer(int N,int kk)
{
  double result;

  result=2.0/((double)N*N)*(1-fabs(2.0*kk+1)/((double)N));

  return result;
}

/** Computes the damping factor for the modified Jackson kernel.
 */
double modified_jackson2(int N,int kk)
{
  int kj;
  double result;
  double n=(N/2+1)/2;
  double k;
  
  for(result=0,kj=kk;kj<=kk+1;kj++)
    {
      k=fabs(kj);
      if(k/n<1)
	result+= 1 - (3.0*k + 6.0*n*pow(k,2) - 3.0*pow(k,3)) 
	  / (2.0*n*(2.0*pow(n,2)+1.0));
      else
	result+= (2*n-k)*(pow(2*n-k,2)-1) / (2.0*n*(2.0*pow(n,2)+1.0));
    }
      
  return result;
}

/** Computes the damping factor for the modified generalised Jackson kernel.
 */
double modified_jackson4(int N,int kk)
{
  int kj;
  double result;
  double n=(N/2+3)/4;
  double k;
  double normalisation=(2416*pow(n,7)+1120*pow(n,5)+784*pow(n,3)+720*n);
  
  for(result=0,kj=kk;kj<=kk+1;kj++)
    {
      k=fabs(kj);

      if(k/n<1)
	result+= 1 - (1260*k + (1680*pow(n,5)+2240*pow(n,3)+2940*n)*pow(k,2) -
		      1715*pow(k,3) - (560*pow(n,3)+1400*n)*pow(k,4) + 490*
		      pow(k,5) + 140*n*pow(k,6) - 35*pow(k,7)) / normalisation;
      
      if((1<=k/n)&&(k/n<2))
  	result+= ((2472*pow(n,7)+336*pow(n,5)+3528*pow(n,3)-1296*n) - 
		  (392*pow(n,6)-3920*pow(n,4)+8232*pow(n,2)-756)*k -
		  (504*pow(n,5)+10080*pow(n,3)-5292*n)*pow(k,2) - 
		  (1960*pow(n,4)-7840*pow(n,2)+1029)*pow(k,3) +
		  (2520*pow(n,3)-2520*n)*pow(k,4) - (1176*pow(n,2)-294)*pow(k,5)
		  + 252*n*pow(k,6) - 21*pow(k,7)) / normalisation;
      
      if((2<=k/n)&&(k/n<3))
 	result+= (-(1112*pow(n,7)-12880*pow(n,5)+7448*pow(n,3)-720*n) +
		  (12152*pow(n,6)-27440*pow(n,4)+8232*pow(n,2)-252)*k -
		  (19320*pow(n,5)-21280*pow(n,3)+2940*n)*pow(k,2) +
		  (13720*pow(n,4)-7840*pow(n,2)+343)*pow(k,3) -
		  (5320*pow(n,3)-1400*n)*pow(k,4) + (1176*pow(n,2)-98)*pow(k,5)
		  - 140*n*pow(k,6) + 7*pow(k,7)) / normalisation;
      
      if((3<=k/n)&&(k/n<4))
	result+= ((4*n-k)*(pow(4*n-k,2)-1)*(pow(4*n-k,2)-4)*(pow(4*n-k,2)-9)) /
	  normalisation;
    }

  return result;
}

/** Computes the damping factor for the modified Sobolev kernel.
 */
double modified_sobolev(double mu,int kk)
{
  double result;
  int kj,k;
  
  for(result=0,kj=kk;kj<=kk+1;kj++)
    {
      k=fabs(kj);
      if(k==0)
	result+=1;
      else
	result+=pow(k,-2*mu);
    }

  return result;
}

/** Computes the damping factor for the modified multiquadric kernel.
 */
double modified_multiquadric(double mu,double c,int kk)
{
  double result;
  int kj,k;
  
  for(result=0,kj=kk;kj<=kk+1;kj++)
    {
      k=fabs(kj);
      result+=pow(k*k+c*c,-mu);
    }

  return result;
}
