/* window functions -----------------------------------------------------------*/
#ifdef DIRAC_DELTA
#define PHI_HUT(k,d) 1.0
#define PHI(x,d) (fabs((x))<10e-8)? 1.0 : 0.0
#define WINDOW_HELP_INIT(d)
#define WINDOW_HELP_FINALIZE
#define WINDOW_HELP_ESTIMATE_m {ths->m = 0;}
#endif

#ifdef GAUSSIAN
#define PHI_HUT(k,d) ((double)exp(-(pow(PI*(k)/ths->n[d],2.0)*           \
                               ths->b[d])))
#define PHI(x,d) ((double)exp(-pow((x)*ths->n[d],2.0)/                   \
                           ths->b[d])/sqrt(PI*ths->b[d]))
#define WINDOW_HELP_INIT {                                                     \
 int idx;                                                                      \
 ths->b = (double*) fftw_malloc(ths->d*sizeof(double));            \
 for(idx=0; idx<ths->d; idx++)                                           \
  ths->b[idx]=((double)2*ths->sigma[idx])/                         \
   (2*ths->sigma[idx]-1)*(((double)ths->m) / PI);                  \
} 
#define WINDOW_HELP_FINALIZE {fftw_free(ths->b);}
#define WINDOW_HELP_ESTIMATE_m {ths->m =12;}
#endif

#ifdef KAISER_BESSEL
#define PHI_HUT(k,d) ((double)nfft_i0( ths->m*sqrt( pow(ths->b[d],2) -  \
                               pow(2*PI*(k)/ths->n[d],2)))) 
#define PHI(x,d) ((double)((pow(ths->m,2)-pow((x)*ths->n[d],2))>0)?\
                   sinh(ths->b[d]*sqrt(pow(ths->m,2)-              \
                   pow((x)*ths->n[d],2)))/(PI*sqrt(pow(ths->m,2)-  \
                   pow((x)*ths->n[d],2))): (((pow(ths->m,2)-       \
                   pow((x)*ths->n[d],2))<0)? sin(ths->b[d]*        \
                   sqrt(pow(ths->n[d]*(x),2)-pow(ths->m,2)))/      \
                   (PI*sqrt(pow(ths->n[d]*(x),2)-pow(ths->m,2))):  \
                   1.0))
#define WINDOW_HELP_INIT {                                                     \
 int idx;                                                                      \
 ths->b = (double*) fftw_malloc(ths->d*sizeof(double));            \
 for(idx=0; idx<ths->d; idx++)                                           \
  ths->b[idx] = ((double)PI*(2.0-1.0/ths->sigma[idx]));            \
}
#define WINDOW_HELP_FINALIZE {fftw_free(ths->b);}
#define WINDOW_HELP_ESTIMATE_m {ths->m = 6;}
#endif

#ifdef B_SPLINE
#define PHI_HUT(k,d) ((double)(((k)==0)? 1.0/ths->n[(d)] :                 \
                       pow(sin((k)*PI/ths->n[(d)])/((k)*PI/ths->n[(d)])\
                       ,2*ths->m)/ths->n[(d)]))
#define PHI(x,d) (nfft_bspline(2*ths->m,((x)*ths->n[(d)])+                \
                 (double)ths->m,ths->spline_coeffs)/ths->n[(d)])
/* wo die /n herkommt ??????????????????? */
#define WINDOW_HELP_INIT {                                                     \
 ths->spline_coeffs= (double*)fftw_malloc(2*ths->m*sizeof(double));\
}
#define WINDOW_HELP_FINALIZE {fftw_free(ths->spline_coeffs);}
#define WINDOW_HELP_ESTIMATE_m {ths->m =11;}
#endif

#ifdef SINC_POWER
#define PHI_HUT(k,d) (nfft_bspline(2*ths->m,((double)2*ths->m*(k))/     \
                      ((2*ths->sigma[(d)]-1)*ths->n[(d)]/ths->sigma[(d)])+             \
                      (double)ths->m,ths->spline_coeffs))
#define PHI(x,d) ((double)(ths->n[(d)]/ths->sigma[(d)]*(2*ths->sigma[(d)]-1)/          \
                  (2*ths->m)*pow(nfft_sinc(PI*ths->n[(d)]/ths->sigma[(d)]*(x)*            \
                  (2*ths->sigma[(d)]-1)/(2*ths->m)),                 \
                  2*ths->m)/ths->n[(d)]))
#define WINDOW_HELP_INIT {                                                     \
 ths->spline_coeffs= (double*)fftw_malloc(2*ths->m*sizeof(double));\
}
#define WINDOW_HELP_FINALIZE {fftw_free(ths->spline_coeffs);}
#define WINDOW_HELP_ESTIMATE_m {ths->m = 9;}
#endif
