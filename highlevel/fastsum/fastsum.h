#ifndef fastsum_h_inc
#define fastsum_h_inc

#include <complex.h>
#include "nfft3.h"


/* OPEN */
Zur schnellen Summation:
0) Bekommt Init einen Funktionspointer auf die Kernfunktion?
   (Der Grund für die Defines bei der NFFT ist, dass PHI innerhalb einer
    Schleife während der eigentlichen Berechnungsphase benötigt wird.
    Hier wird der Kern aber hauptsächlich zur Initialisierung benötigt, oder?)
1) verschiedene N_0, N_1, ... sind nur sinnvoll für K(x)=\phi(x^T \Sigma x),
   d.h. nicht für radiale, oder?
2) die Knoten lagen in [-1/4+eps/2,1/4-eps/2], da am Rand regularisiert wurde?
3) das Lösen von Gleichungssystemen K(x_j-y_l) soll mit vorgesehen werden?
   Dann Restriktionen für Namen und adjoint notwendig.
4) beim Gauss haben wir die Periode p?
5) die 2Pkt-Taylorsumme ist die Regul. der Wahl?
6) Nahfeld - Datenstruktur?
7) fastsum_trafo berechnet auch die Nahfeldbeiträge oder nicht?
8) Sind die Kerne symmetrisch

/** 
 * Constant symbols for precomputation and memory usage (inverse problem)
 */
#define DS_PRE       (1U<< 0) /* precompute all entries of K(x_j-y_l) */
#define FS_NDFT      (1U<< 1) /* approximate K by K_n+Nearfield but multiply with k_n(x_j-y_l) by NDFT */
#define FS_SINGULAR  (1U<< 2) /* regularise singularity at the origin, set up nearfield */
#define FS_BOUNDARY  (1U<< 3) /* regularise boundary, restriction on nodes x_j,y_l */

typedef struct fastsum_plan_ 
{ 
/** api */

  int N_total;                          /**< number of source knots          */
  int M_total;                          /**< number of target knots          */

  complex *alpha;                       /**< source coefficients             */
  complex *f;                           /**< target evaluations              */

  double *x;                            /**< source knots in [-1/4,1/4]      */
  double *y;                            /**< target knots in [-1/4,1/4]      */
 
  unsigned flags;                       /**< flags precomp. and approx.type  */

/** internal */

/** DS_PRE - direct summation */ 
  complex *pre_K;                       /**< precomputed K(x_j-y_l)          */

/** FS__ - fast summation */
  int n;                                /**< expansion degree                */
  complex *b;                           /**< expansion coefficients          */

  double p;                             /**< period, used in gauss trafo???? */

  nfft_plan *mv1;                       /**< source nfft plan                */
  nfft_plan *mv2;                       /**< target nfft plan                */

/* near field??? */

/* things for computing *b - are they used only once?? */
  fftw_plan *fft;
} fastsum_plan;

/** initialisation for the fast summation, assumes regular kernel
 */
void fastsum_init_simple(fastsum_plan *ths, complex (*kernel)());

/** initialisation for the fast summation
 */
void fastsum_init_advanced(fastsum_plan *ths, complex (*kernel)(), unsigned flags);

/** initialisation for the fast summation
 */
void fastsum_init_guru(fastsum_plan *ths, nfft_plan *mv, complex (*kernel)(), unsigned flags);

/** node dependent precomputation step, i.e., calls nfft_precompute_one_psi()
 */
void fastsum_precompute(fastsum_plan *ths);

/** computes fast summation
 */
void fastsum_trafo(fastsum_plan *ths);

/** finalisation for the fast summation
 */
void fastsum_finalize(fastsum_plan *ths);
/** @} 
 */ 

#endif
