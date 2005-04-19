/* window functions -----------------------------------------------------------*/
#ifdef DIRAC_DELTA
#define PHI_HUT(k,d) 1.0
#define PHI(x,d) (fabs((x))<10e-8)? 1.0 : 0.0
#define WINDOW_HELP_INIT(d)
#define WINDOW_HELP_FINALIZE
#define WINDOW_HELP_ESTIMATE_m {this->m = 0;}
#endif

#ifdef GAUSSIAN
#define PHI_HUT(k,d) ((double)exp(-(pow(PI*(k)/this->n[d],2.0)*           \
                               this->b[d])))
#define PHI(x,d) ((double)exp(-pow((x)*this->n[d],2.0)/                   \
                           this->b[d])/sqrt(PI*this->b[d]))
#define WINDOW_HELP_INIT {                                                     \
 int idx;                                                                      \
 this->b = (double*) fftw_malloc(this->d*sizeof(double));            \
 for(idx=0; idx<this->d; idx++)                                           \
  this->b[idx]=((double)2*this->sigma[idx])/                         \
   (2*this->sigma[idx]-1)*(((double)this->m) / PI);                  \
} 
#define WINDOW_HELP_FINALIZE {fftw_free(this->b);}
#define WINDOW_HELP_ESTIMATE_m {this->m =12;}
#endif

#ifdef KAISER_BESSEL
#define PHI_HUT(k,d) ((double)i0( this->m*sqrt( pow(this->b[d],2) -  \
                               pow(2*PI*(k)/this->n[d],2)))) 
#define PHI(x,d) ((double)((pow(this->m,2)-pow((x)*this->n[d],2))>0)?\
                   sinh(this->b[d]*sqrt(pow(this->m,2)-              \
                   pow((x)*this->n[d],2)))/(PI*sqrt(pow(this->m,2)-  \
                   pow((x)*this->n[d],2))): (((pow(this->m,2)-       \
                   pow((x)*this->n[d],2))<0)? sin(this->b[d]*        \
                   sqrt(pow(this->n[d]*(x),2)-pow(this->m,2)))/      \
                   (PI*sqrt(pow(this->n[d]*(x),2)-pow(this->m,2))):  \
                   1.0))
#define WINDOW_HELP_INIT {                                                     \
 int idx;                                                                      \
 this->b = (double*) fftw_malloc(this->d*sizeof(double));            \
 for(idx=0; idx<this->d; idx++)                                           \
  this->b[idx] = ((double)PI*(2.0-1.0/this->sigma[idx]));            \
}
#define WINDOW_HELP_FINALIZE {fftw_free(this->b);}
#define WINDOW_HELP_ESTIMATE_m {this->m = 6;}
#endif

#ifdef B_SPLINE
#define PHI_HUT(k,d) ((double)(((k)==0)? 1.0/this->n[(d)] :                 \
                       pow(sin((k)*PI/this->n[(d)])/((k)*PI/this->n[(d)])\
                       ,2*this->m)/this->n[(d)]))
#define PHI(x,d) (bspline(2*this->m,((x)*this->n[(d)])+                \
                 (double)this->m,this->spline_coeffs)/this->n[(d)])
/* wo die /n herkommt ??????????????????? */
#define WINDOW_HELP_INIT {                                                     \
 this->spline_coeffs= (double*)fftw_malloc(2*this->m*sizeof(double));\
}
#define WINDOW_HELP_FINALIZE {fftw_free(this->spline_coeffs);}
#define WINDOW_HELP_ESTIMATE_m {this->m =11;}
#endif

#ifdef SINC_POWER
#define PHI_HUT(k,d) (bspline(2*this->m,((double)2*this->m*(k))/     \
                      ((2*this->sigma[(d)]-1)*this->N[(d)])+             \
                      (double)this->m,this->spline_coeffs))
#define PHI(x,d) ((double)(this->N[(d)]*(2*this->sigma[(d)]-1)/          \
                  (2*this->m)*pow(sinc(PI*this->N[(d)]*(x)*            \
                  (2*this->sigma[(d)]-1)/(2*this->m)),                 \
                  2*this->m)/this->n[(d)]))
#define WINDOW_HELP_INIT {                                                     \
 this->spline_coeffs= (double*)fftw_malloc(2*this->m*sizeof(double));\
}
#define WINDOW_HELP_FINALIZE {fftw_free(this->spline_coeffs);}
#define WINDOW_HELP_ESTIMATE_m {this->m = 9;}
#endif
