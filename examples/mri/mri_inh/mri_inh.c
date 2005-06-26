//#include "solver.c"
#include "mri_inh.h"
#include "options.h"
#include "util.h"
#include "window_defines.h"
#include "math.h"

#define PHI_periodic(x) ((x>0.5)?(PHI(x-1.0,0)):((x<-0.5)?PHI(x+1.0,0):PHI(x,0)))

void mri_inh_trafo(mri_inh_plan *that) {
  int l,j;
  complex *f = (complex*) fftw_malloc(that->M_total*sizeof(complex));
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_1d(ths,that->N3,1);
  memset(f,0,that->M_total*sizeof(complex));
  for(j=0;j<that->N_total;j++)
  {
    that->f_hat[j]/=PHI_HUT(ths->n[0]*that->w[j],0);
    f_hat[j]=that->f_hat[j];
  }

  for(l=-ths->n[0]/2;l<ths->n[0]/2;l++) {
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]*=cexp(-2*PI*I*that->w[j]*l);
    nfft_trafo((nfft_plan*)that);
    for(j=0;j<that->M_total;j++)
      f[j]+=that->f[j]*PHI_periodic(that->t[j]-((double)l)/((double)ths->n[0]));
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]=f_hat[j];
  }

  fftw_free(that->f);
  that->f=f;
  fftw_free(f_hat);
  
  nfft_finalize(ths);
  free(ths);
}

void mri_inh_adjoint(mri_inh_plan *that) {
  int l,j;
  complex *f = (complex*) fftw_malloc(that->M_total*sizeof(complex));
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_1d(ths,that->N3,1);
  memset(f_hat,0,that->N_total*sizeof(complex));
  for(j=0;j<that->M_total;j++)
  {
    f[j]=that->f[j];
  }


  
  for(l=-ths->n[0]/2;l<ths->n[0]/2;l++) {
    
    for(j=0;j<that->M_total;j++)
      that->f[j]*=PHI_periodic(that->t[j]-((double)l)/((double)ths->n[0]));
    printf("%i -\n",l);
    nfft_adjoint((nfft_plan*)that);
    printf("%i +\n",l);
    for(j=0;j<that->N_total;j++)
      f_hat[j]+=that->f_hat[j]*cexp(2*PI*I*that->w[j]*l);
    for(j=0;j<that->M_total;j++)
      that->f[j]=f[j];
  }

  for(j=0;j<that->N_total;j++)
  {
    f_hat[j] /= PHI_HUT(ths->n[0]*that->w[j],0);
  }



  fftw_free(that->f_hat);
  that->f_hat=f_hat;
  fftw_free(f);
  
  nfft_finalize(ths);
  free(ths);
}

void mri_inh_init_guru(mri_inh_plan *ths, int d, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags) {

  nfft_init_guru((nfft_plan*)ths,d,N,M,n,m,nfft_flags,fftw_flags);
  ths->N3=N[2];
  ths->t = (double*) fftw_malloc(ths->M_total*sizeof(double));
  ths->w = (double*) fftw_malloc(ths->N_total*sizeof(double));
}

void mri_inh_finalize(mri_inh_plan *ths) {
  //fftw_free(ths->t);
  fftw_free(ths->w);
  nfft_finalize((nfft_plan*)ths);
}


void mri_inh_trafo2(mri_inh_plan *that) {

}

void mri_inh_adjoint2(mri_inh_plan *that) {
  int l,j;
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_1d(ths,that->N3,1);
  memset(f_hat,0,that->N_total*sizeof(complex));

  for(j=0;j<that->M_total;j++)
  {
    that->f[j] /= PHI_HUT(ths->n[0]*that->x[3*j+2],0);
  }
  
  nfft_adjoint((nfft_plan*)that);
  

  for(j=0;j<that->N[0]*that->N[1];j++) {
    for(l=-ths->n[0]/2;l<ths->n[0]/2;l++)
    {
      f_hat[j]+= that->f_hat[j*ths->n[0]+(l+ths->n[0]/2)]*PHI_periodic(that->w[j]-((double)l)/((double)ths->n[0]));
    }
  }





  fftw_free(that->f_hat);
  that->f_hat=f_hat;
  
  nfft_finalize(ths);
  free(ths);
}

void mri_inh_init_guru2(mri_inh_plan *ths, int d, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags) {
  ths->N3=N[2];
  N[2]=n[2];
  nfft_init_guru((nfft_plan*)ths,d,N,M,n,m,nfft_flags,fftw_flags);
  ths->w = (double*) fftw_malloc(ths->N_total*sizeof(double));
}

void mri_inh_finalize2(mri_inh_plan *ths) {
  fftw_free(ths->w);
  nfft_finalize((nfft_plan*)ths);
}


#define MACRO_SOLVER_IMPL(MV, FLT)                                            \
                                                                              \
F(MV, FLT, init_advanced, i ## MV ## _plan *ths, MV ## _plan *mv,             \
                          unsigned i ## MV ## _flags)                         \
{                                                                             \
  ths->mv = mv;                                                               \
  ths->flags = i ## MV ## _flags;                                             \
                                                                              \
  ths->y          = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));          \
  ths->r_iter     = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));          \
  ths->f_hat_iter = (FLT*)fftw_malloc(ths->mv->N_total*sizeof(FLT));          \
  ths->p_hat_iter = (FLT*)fftw_malloc(ths->mv->N_total*sizeof(FLT));          \
                                                                              \
  if(ths->flags & LANDWEBER)                                                  \
    ths->z_hat_iter = ths->p_hat_iter;                                        \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    {                                                                         \
      ths->z_hat_iter = ths->p_hat_iter;                                      \
      ths->v_iter     = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));      \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    {                                                                         \
      ths->z_hat_iter = (FLT*)fftw_malloc(ths->mv->N_total*sizeof(FLT));      \
      ths->v_iter     = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));      \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    ths->z_hat_iter = ths->p_hat_iter;                                        \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->w = (double*) fftw_malloc(ths->mv->M_total*sizeof(double));          \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->w_hat = (double*) fftw_malloc(ths->mv->N_total*sizeof(double));      \
}                                                                             \
                                                                              \
/** void i<mv>_init */                                                        \
F(MV, FLT, init, i ## MV ## _plan *ths, MV ## _plan *mv)                      \
{                                                                             \
  i ## MV ## _init_advanced(ths, mv, CGNR);                                   \
} /* void i<mv>_init */                                                       \
                                                                              \
/** void i<mv>_before_loop */                                                 \
F(MV, FLT, before_loop,   i ## MV ## _plan *ths)                              \
{                                                                             \
  cp_ ## FLT(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);              \
                                                                              \
  SWAP_ ## FLT(ths->r_iter, ths->mv->f);                                      \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->r_iter, ths->mv->f);                                      \
                                                                              \
  upd_axpy_ ## FLT(ths->r_iter, -1.0, ths->y, ths->mv->M_total);              \
                                                                              \
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))       \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_WEIGHT)                                      \
        ths->dot_r_iter =                                                     \
            dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);             \
      else                                                                    \
        ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);         \
    }                                                                         \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                              \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                              \
                                                                              \
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))       \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_DAMP)                                        \
        ths->dot_z_hat_iter =                                                 \
          dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);       \
      else                                                                    \
        ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total); \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;                                \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    cp_ ## FLT(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);           \
} /* void i<mv>_before_loop */                                                \
                                                                              \
/** void i<mv>_loop_one_step_landweber */                                     \
F(MV, FLT, loop_one_step_landweber, i ## MV ## _plan *ths)                    \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
                      ths->z_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,       \
                     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  cp_ ## FLT(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);              \
                                                                              \
  SWAP_ ## FLT(ths->r_iter,ths->mv->f);                                       \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->r_iter,ths->mv->f);                                       \
                                                                              \
  upd_axpy_ ## FLT(ths->r_iter, -1.0, ths->y, ths->mv->M_total);              \
                                                                              \
  if(ths->flags & NORMS_FOR_LANDWEBER)                                        \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_WEIGHT)                                      \
        ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter,ths->w, ths->mv->M_total);\
      else                                                                    \
        ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);         \
    }                                                                         \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
                                                                              \
  if(ths->flags & NORMS_FOR_LANDWEBER)                                        \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_DAMP)                                        \
        ths->dot_z_hat_iter =                                                 \
          dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);       \
      else                                                                    \
        ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total); \
    }                                                                         \
} /* void i<mv>_loop_one_step_landweber */                                    \
                                                                              \
/** void i<mv>_loop_one_step_steepest_descent */                              \
F(MV, FLT, loop_one_step_steepest_descent, i ## MV ## _plan *ths)             \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter,                 \
                 ths->mv->N_total);                                           \
  else                                                                        \
    cp_ ## FLT(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);            \
                                                                              \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_v_iter = dot_w_ ## FLT(ths->v_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_v_iter = dot_ ## FLT(ths->v_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
                      ths->z_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,       \
                     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->v_iter,                \
                   ths->mv->M_total);                                         \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_z_hat_iter =                                                     \
      dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);           \
  else                                                                        \
    ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);     \
} /* void i<mv>_loop_one_step_steepest_descent */                             \
                                                                              \
/** void i<mv>_loop_one_step_cgnr */                                          \
F(MV, FLT, loop_one_step_cgnr, i ## MV ## _plan *ths)                         \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,                 \
                 ths->mv->N_total);                                           \
  else                                                                        \
    cp_ ## FLT(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);            \
                                                                              \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_v_iter = dot_w_ ## FLT(ths->v_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_v_iter = dot_ ## FLT(ths->v_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
                      ths->p_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,       \
                     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->v_iter,                \
                   ths->mv->M_total);                                         \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
                                                                              \
  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_z_hat_iter =                                                     \
      dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);           \
  else                                                                        \
    ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);     \
                                                                              \
  /*-----------------*/                                                       \
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;             \
                                                                              \
  /*-----------------*/                                                       \
  upd_axpy_ ## FLT(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,          \
                   ths->mv->N_total);                                         \
} /* void i<mv>_loop_one_step_cgnr */                                         \
                                                                              \
/** void i<mv>_loop_one_step_cgne */                                          \
F(MV, FLT, loop_one_step_cgne, i ## MV ## _plan *ths)                         \
{                                                                             \
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
                      ths->p_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,       \
                     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,                 \
                 ths->mv->N_total);                                           \
  else                                                                        \
    cp_ ## FLT(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);            \
                                                                              \
  MV ## _trafo(ths->mv);                                                      \
                                                                              \
  upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->mv->f,                 \
                   ths->mv->M_total);                                         \
                                                                              \
  ths->dot_r_iter_old = ths->dot_r_iter;                                      \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;                     \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  MV ## _adjoint(ths->mv);                                                    \
                                                                              \
  upd_axpy_ ## FLT(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,           \
                   ths->mv->N_total);                                         \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_p_hat_iter =                                                     \
      dot_w_ ## FLT(ths->p_hat_iter, ths->w_hat, ths->mv->N_total);           \
  else                                                                        \
    ths->dot_p_hat_iter = dot_ ## FLT(ths->p_hat_iter, ths->mv->N_total);     \
} /* void i<mv>_loop_one_step_cgne */                                         \
                                                                              \
/** void i<mv>_loop_one_step */                                               \
F(MV, FLT, loop_one_step, i ## MV ## _plan *ths)                              \
{                                                                             \
  if(ths->flags & LANDWEBER)                                                  \
    i ## MV ## _loop_one_step_landweber(ths);                                 \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    i ## MV ## _loop_one_step_steepest_descent(ths);                          \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    i ## MV ## _loop_one_step_cgnr(ths);                                      \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    i ## MV ## _loop_one_step_cgne(ths);                                      \
} /* void i<mv>_loop_one_step */                                              \
                                                                              \
/** void i<mv>_finalize */                                                    \
F(MV, FLT, finalize, i ## MV ## _plan *ths)                                   \
{                                                                             \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    fftw_free(ths->w);                                                        \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    fftw_free(ths->w_hat);                                                    \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    {                                                                         \
      fftw_free(ths->v_iter);                                                 \
      fftw_free(ths->z_hat_iter);                                             \
    }                                                                         \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    fftw_free(ths->v_iter);                                                   \
                                                                              \
  fftw_free(ths->p_hat_iter);                                                 \
  fftw_free(ths->f_hat_iter);                                                 \
                                                                              \
  fftw_free(ths->r_iter);                                                     \
  fftw_free(ths->y);                                                          \
} /** void i<mv>_finalize */

MACRO_SOLVER_IMPL(mri_inh, complex)