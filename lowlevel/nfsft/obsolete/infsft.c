#include "infsft.h"
#include "api.h"
#include "util.h"
#include <stdlib.h>

void infsft_init_old(infsft_plan_old this_iplan, nfsft_plan_old direct_plan)
{
  infsft_init_guru_old(this_iplan,direct_plan,NFSFT_CGNR_E);
}

infsft_plan_old infsft_make_plan_old()
{
  return (infsft_plan_old) malloc(sizeof(struct infsft_plan_s_old));
}

void infsft_init_guru_old(infsft_plan_old this_iplan, nfsft_plan_old direct_plan, 
                      int infsft_flags_old)
{
  int n;
  
  this_iplan->direct_plan = direct_plan;
  this_iplan->infsft_flags_old = infsft_flags_old;
  
  this_iplan->given_f = (complex*)
    fftw_malloc(this_iplan->direct_plan->D*sizeof(complex));
  
  this_iplan->r_iter = (complex*)
    fftw_malloc(this_iplan->direct_plan->D*sizeof(complex));
  
  this_iplan->f_hat_iter = (complex**)
    fftw_malloc((2*this_iplan->direct_plan->M+1)*sizeof(complex*));
  for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
  {
    this_iplan->f_hat_iter[n] = (complex*)
    fftw_malloc((this_iplan->direct_plan->N+1)*sizeof(complex));    
  }
  
  this_iplan->p_hat_iter = (complex**)
    fftw_malloc((2*this_iplan->direct_plan->M+1)*sizeof(complex*));
  for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
  {
    this_iplan->p_hat_iter[n] = (complex*)
    fftw_malloc((this_iplan->direct_plan->N+1)*sizeof(complex));    
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_LANDWEBER)
  {
    this_iplan->z_hat_iter = this_iplan->p_hat_iter;
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_STEEPEST_DESCENT)
  {
    this_iplan->z_hat_iter = this_iplan->p_hat_iter;
    
    this_iplan->v_iter = (complex*)
      fftw_malloc(this_iplan->direct_plan->D*sizeof(complex));
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNR_E)
  {
    this_iplan->z_hat_iter = (complex**)
      fftw_malloc((2*this_iplan->direct_plan->M+1)*sizeof(complex*));
    for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
    {
      this_iplan->z_hat_iter[n] = (complex*)
      fftw_malloc((this_iplan->direct_plan->N+1)*sizeof(complex));    
    }
    
    this_iplan->v_iter = (complex*)
      fftw_malloc(this_iplan->direct_plan->D*sizeof(complex));
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNE_R)
  {
    this_iplan->z_hat_iter = this_iplan->p_hat_iter;
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_ITERATE_2nd)
  {
    this_iplan->f_hat_iter_2nd = (complex**) 
      fftw_malloc((2*this_iplan->direct_plan->M+1)*sizeof(fftw_complex*));
    for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
    {
      this_iplan->f_hat_iter_2nd[n] = (complex*)
      fftw_malloc((this_iplan->direct_plan->N+1)*sizeof(complex));    
    }
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
  {
    this_iplan->w = 
      (double*) fftw_malloc(this_iplan->direct_plan->D*sizeof(double));
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
  {
    this_iplan->w_hat = 
      (double**) fftw_malloc((2*this_iplan->direct_plan->M+1)*sizeof(double*));
    for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
    {
      this_iplan->w_hat[n] = (double*)
      fftw_malloc((this_iplan->direct_plan->N+1)*sizeof(double));    
    }
  }

  //----
  /*copyc_hat_old(this_iplan->f_hat_iter,this_iplan->direct_plan->f_hat,
            this_iplan->direct_plan->M);
  copyc(this_iplan->given_f,this_iplan->direct_plan->f,
            this_iplan->direct_plan->D);*/
  //----
}

void infsft_before_loop_help_old(infsft_plan_old this_iplan)
{
  /** step 2
   *  overwrites this_iplan->direct_plan->f_hat 
   *  overwrites this_iplan->r_iter
   */
  
  complex *temp;
  complex **temp_hat;
  
  //vpr_c_hat(this_iplan->direct_plan->f_hat, this_iplan->direct_plan->M, "f_hat");
  //vpr_c_hat(this_iplan->f_hat_iter, this_iplan->direct_plan->M, "f_hat_iter");

  copyc_hat_old(this_iplan->direct_plan->f_hat, this_iplan->f_hat_iter,
        this_iplan->direct_plan->M);
  
  //vpr_c_hat(this_iplan->direct_plan->f_hat, this_iplan->direct_plan->M, "f_hat");
  //vpr_c_hat(this_iplan->f_hat_iter, this_iplan->direct_plan->M, "f_hat_iter");
    
  SWAPCT(this_iplan->r_iter,this_iplan->direct_plan->f,temp)
          
  nfsft_trafo_old(this_iplan->direct_plan);

  SWAPCT(this_iplan->r_iter,this_iplan->direct_plan->f,temp)
    
  /*printf("D = %d\n",this_iplan->direct_plan->D);
  vpr_c(this_iplan->r_iter, this_iplan->direct_plan->D, "r_iter");
  vpr_c(this_iplan->given_f, this_iplan->direct_plan->D, "given_f");*/

  updatec_axpy_2_old(this_iplan->r_iter, -1.0, this_iplan->given_f,
               this_iplan->direct_plan->D);
  
  printf("D = %d\n",this_iplan->direct_plan->D);
  //vpr_c(this_iplan->r_iter, this_iplan->direct_plan->D, "r_iter");
  //vpr_c(this_iplan->given_f, this_iplan->direct_plan->D, "given_f");
  
  if((!(this_iplan->infsft_flags_old & NFSFT_LANDWEBER)) ||
     (this_iplan->infsft_flags_old & NFSFT_NORMS_FOR_LANDWEBER))
  {
    if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
    {  
      //fprintf(stderr,"with weight!\n");
      this_iplan->dot_r_iter =
        dot_w_complex(this_iplan->r_iter, this_iplan->w,
                      this_iplan->direct_plan->D);
    }
    else
    {
      //fprintf(stderr,"without weight!\n");
      this_iplan->dot_r_iter =
        dot_complex(this_iplan->r_iter, this_iplan->direct_plan->D);
    }
  }
  else
  {
    //fprintf(stderr,"no dot_r_iter!\n");
  }
  //fflush(stderr);
  
  //fprintf(stderr,"dot_r_iter = %f\n",this_iplan->dot_r_iter);
  
  /** step 3
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->z_hat_iter resp. this_iplan->z_hat_iter
    */ 
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
    cp_w_complex(this_iplan->direct_plan->f, this_iplan->w, this_iplan->r_iter,
            this_iplan->direct_plan->D);
  else
    cp_complex(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->D);
    
  SWAPCT(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat, temp_hat);
  nfsft_adjoint_old(this_iplan->direct_plan);
  //vpr_c(this_iplan->direct_plan->f,this_iplan->direct_plan->D,"f");
  //vpr_c_hat(this_iplan->direct_plan->f_hat,this_iplan->direct_plan->M,"f_hat");
  SWAPCT(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat, temp_hat);
   
  if((!(this_iplan->infsft_flags_old & NFSFT_LANDWEBER)) ||
     (this_iplan->infsft_flags_old & NFSFT_NORMS_FOR_LANDWEBER))
  {
    if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
      this_iplan->dot_z_hat_iter = 
        dotproductc_w_hat_old(this_iplan->z_hat_iter, this_iplan->w_hat,
                          this_iplan->direct_plan->M);
    else
    {
      this_iplan->dot_z_hat_iter =
        dotproductc_hat_old(this_iplan->z_hat_iter, this_iplan->direct_plan->M);
      //fprintf(stderr,"dot_z_hat_iter\n");
    }
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNE_R)
    this_iplan->dot_p_hat_iter = this_iplan->dot_z_hat_iter;
} /* void infft_before_loop_help */

void infsft_before_loop_old(infsft_plan_old this_iplan)
{
  infsft_before_loop_help_old(this_iplan);
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNR_E)
  {
    /** step 4-6
    *  overwrites this_iplan->f_hat_iter_2nd
    */
    if(this_iplan->infsft_flags_old & NFSFT_ITERATE_2nd)
      copyc_hat_old(this_iplan->f_hat_iter_2nd, this_iplan->f_hat_iter,
            this_iplan->direct_plan->M);
    
    /** step 7
    *  overwrites this_iplan->p_hat_iter
    */
    copyc_hat_old(this_iplan->p_hat_iter, this_iplan->z_hat_iter,
          this_iplan->direct_plan->M);
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNE_R)
  {
    /** step 4-7
    *  overwrites this_iplan->f_hat_iter_2nd
    */
    if(this_iplan->infsft_flags_old & NFSFT_ITERATE_2nd)
    {
      this_iplan->gamma_iter=1.0;
      copyc_hat_old(this_iplan->f_hat_iter_2nd, this_iplan->f_hat_iter,
            this_iplan->direct_plan->M);
    }
  }
} /* void infft_before_loop */

/*void infft_loop_one_step_landweber(infft_plan *this_iplan)
{*/
  /** step 5
  *  updates this_iplan->f_hat_iter
  */
  /*if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                  this_iplan->w_hat, this_iplan->z_hat_iter,
                  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                 this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);*/
  
  /** step 6
    *  original residual, not the updated residual,
    *  overwrites this_iplan->r_iter
    *  overwrites this_iplan->direct_plan->f_hat 
    */
  /*copyc(this_iplan->direct_plan->f_hat, this_iplan->f_hat_iter,
        this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->r_iter,this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->r_iter,this_iplan->direct_plan->f);
  
  updatec_axpy(this_iplan->r_iter, -1.0, this_iplan->given_f, 
               this_iplan->direct_plan->M);
  
  if(this_iplan->infft_flags & NFSFT_NORMS_FOR_LANDWEBER)
  {
    if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_WEIGHT)
      this_iplan->dot_r_iter = dotproductc_w(this_iplan->r_iter,this_iplan->w,
                                             this_iplan->direct_plan->M);
    else
      this_iplan->dot_r_iter =
        dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);
  }*/
  
  /** step 7
    *  overwrites this_iplan->direct_plan->f 
    *  overwrites this_iplan->z_hat_iter
    */
  /*if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w,
            this_iplan->r_iter, this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  
  if(this_iplan->infft_flags & NFSFT_NORMS_FOR_LANDWEBER)
  {
    if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_DAMP)
      this_iplan->dot_z_hat_iter = 
        dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat,
                      this_iplan->direct_plan->N_L);
    else
      this_iplan->dot_z_hat_iter =
        dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  }
}*/ /* void infft_loop_one_step_landweber */

/*void infft_loop_one_step_steepest_descent(infft_plan *this_iplan)
{*/
  /** step 5
  *  overwrites this_iplan->direct_plan->f_hat 
  *  overwrites this_iplan->v_iter
  */
  /*if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_DAMP)
    copyc_w(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
            this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  else
    copyc(this_iplan->direct_plan->f_hat, this_iplan->z_hat_iter,
          this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  
  if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_WEIGHT)
    this_iplan->dot_v_iter = dotproductc_w(this_iplan->v_iter, this_iplan->w,
                                           this_iplan->direct_plan->M);
  else
    this_iplan->dot_v_iter =
      dotproductc(this_iplan->v_iter, this_iplan->direct_plan->M);*/
  
  /** step 6
    */
  //this_iplan->alpha_iter =
  //  this_iplan->dot_z_hat_iter / this_iplan->dot_v_iter;
  
  /** step 7
    *  updates this_iplan->f_hat_iter
    */
  /*if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                  this_iplan->w_hat,this_iplan->z_hat_iter,
                  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                 this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);*/
  
  /** step 8
    *  updates this_iplan->r_iter
    */
  /*updatec_xpay(this_iplan->r_iter, -this_iplan->alpha_iter, this_iplan->v_iter,
               this_iplan->direct_plan->M);
  
  if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter = dotproductc_w(this_iplan->r_iter, this_iplan->w,
                                           this_iplan->direct_plan->M);
  else
    this_iplan->dot_r_iter =
      dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);*/
  
  /** step 9
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->z_hat_iter
    */
  /*if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w, this_iplan->r_iter,
            this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  
  if(this_iplan->infft_flags & NFSFT_PRECOMPUTE_DAMP)
    this_iplan->dot_z_hat_iter =
      dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat, 
                    this_iplan->direct_plan->N_L);
  else
    this_iplan->dot_z_hat_iter =
      dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
}*/ /* void infft_loop_one_step_steepest_descent */

void infsft_loop_one_step_cgnr_e(infsft_plan_old this_iplan)
{
  complex *temp;
  complex **temp_hat;
  
  /** step 9
  *  overwrites this_iplan->direct_plan->f_hat 
  *  overwrites this_iplan->v_iter
  */
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
    copyc_w_hat_old(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
            this_iplan->p_hat_iter, this_iplan->direct_plan->M);
  else
    copyc_hat_old(this_iplan->direct_plan->f_hat, this_iplan->p_hat_iter,
          this_iplan->direct_plan->M);
  
  SWAPCT(this_iplan->v_iter,this_iplan->direct_plan->f,temp);
  nfsft_trafo_old(this_iplan->direct_plan);
  SWAPCT(this_iplan->v_iter,this_iplan->direct_plan->f,temp);
  
  //vpr_c(this_iplan->v_iter,this_iplan->direct_plan->D,"v_iter");
  //vpr_c(this_iplan->r_iter,this_iplan->direct_plan->D,"r_iter");
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
  {
    this_iplan->dot_v_iter = dot_w_complex(this_iplan->v_iter, this_iplan->w,
                                           this_iplan->direct_plan->D);
    //fprintf(stderr,"dot_v_iter\n");
  }
  else
    this_iplan->dot_v_iter =
      dot_complex(this_iplan->v_iter, this_iplan->direct_plan->D);
  
  /** step 10
    */
  //vpr_c(this_iplan->v_iter,this_iplan->direct_plan->D,"v_iter");
  //vpr_c_hat(this_iplan->z_hat_iter,this_iplan->direct_plan->M,"z_iter");
  this_iplan->alpha_iter =
    this_iplan->dot_z_hat_iter / this_iplan->dot_v_iter;
  
  //printf("alpha_iter = %f\n",this_iplan->alpha_iter);
  //this_iplan->alpha_iter = 1.0;
  
  /** step 11
    *  updates this_iplan->f_hat_iter
    */
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
    updatec_xpawy_hat_old(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                  this_iplan->w_hat, this_iplan->p_hat_iter,
                      this_iplan->direct_plan->M);
  else
    updatec_xpay_hat_old(this_iplan->f_hat_iter, this_iplan->alpha_iter, 
                 this_iplan->p_hat_iter, this_iplan->direct_plan->M);
  
  /** step 12-15
    *  updates this_iplan->f_hat_iter_2nd
    */
  if (this_iplan->infsft_flags_old & NFSFT_ITERATE_2nd)
  {
    this_iplan->alpha_iter_2nd =
    this_iplan->dot_r_iter / this_iplan->dot_z_hat_iter;
    
    if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
      updatec_xpawy_hat_old(this_iplan->f_hat_iter_2nd, this_iplan->alpha_iter_2nd,
                    this_iplan->w_hat, this_iplan->z_hat_iter,
                    this_iplan->direct_plan->M);
    else
      updatec_xpay_hat_old(this_iplan->f_hat_iter_2nd, this_iplan->alpha_iter_2nd,
                   this_iplan->z_hat_iter, this_iplan->direct_plan->M);
  }
  
  /** step 16
    *  updates this_iplan->r_iter
    */
  //vpr_c(this_iplan->v_iter,this_iplan->direct_plan->D,"v_iter");
  //vpr_c(this_iplan->r_iter,this_iplan->direct_plan->D,"r_iter");

  upd_xpay_complex(this_iplan->r_iter, -this_iplan->alpha_iter, this_iplan->v_iter,
               this_iplan->direct_plan->D);
  
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter = dot_w_complex(this_iplan->r_iter, this_iplan->w,
                                           this_iplan->direct_plan->D);
  else
    this_iplan->dot_r_iter =
      dot_complex(this_iplan->r_iter, this_iplan->direct_plan->D);
  
  /** step 17
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->direct_plan->r_iter
    */
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
    cp_w_complex(this_iplan->direct_plan->f, this_iplan->w, 
            this_iplan->r_iter, this_iplan->direct_plan->D);
  else
    cp_complex(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->D);

  SWAPCT(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat, temp_hat);
  nfsft_adjoint_old(this_iplan->direct_plan);
  SWAPCT(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat, temp_hat);

  this_iplan->dot_z_hat_iter_old = this_iplan->dot_z_hat_iter;
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
    this_iplan->dot_z_hat_iter =
      dotproductc_w_hat_old(this_iplan->z_hat_iter, this_iplan->w_hat, 
                    this_iplan->direct_plan->M);
  else
    this_iplan->dot_z_hat_iter =
      dotproductc_hat_old(this_iplan->z_hat_iter, this_iplan->direct_plan->M);
  
  /** step 18
    */
  this_iplan->beta_iter =
    this_iplan->dot_z_hat_iter / this_iplan->dot_z_hat_iter_old;
  
  /** step 19
    *  updates this_iplan->p_hat_iter
    */
  updatec_axpy_hat_old(this_iplan->p_hat_iter, this_iplan->beta_iter, 
               this_iplan->z_hat_iter, this_iplan->direct_plan->M);
} /* void infft_loop_one_step_cgnr_e */

void infsft_loop_one_step_cgne_r(infsft_plan_old this_iplan)
{
  /** step 9
  */
  this_iplan->alpha_iter =
  this_iplan->dot_r_iter / this_iplan->dot_p_hat_iter;
  
  /** step 10
  *  updates this_iplan->f_hat_iter
  */
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
    updatec_xpawy_hat_old(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                  this_iplan->w_hat, this_iplan->p_hat_iter,
                  this_iplan->direct_plan->M);
  else
    updatec_xpay_hat_old(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                     this_iplan->p_hat_iter, this_iplan->direct_plan->M);
  
  /** step 11
    *  overwrites this_iplan->direct_plan->f_hat 
    *  overwrites this_iplan->direct_plan->f
    *  updates this_iplan->r_iter
    */
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
    copyc_w_hat_old(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
            this_iplan->p_hat_iter, this_iplan->direct_plan->M);
  else
    copyc_hat_old(this_iplan->direct_plan->f_hat, this_iplan->p_hat_iter,
          this_iplan->direct_plan->M);
  
  nfsft_trafo_old(this_iplan->direct_plan);
  
  upd_xpay_complex(this_iplan->r_iter, -this_iplan->alpha_iter,
               this_iplan->direct_plan->f, this_iplan->direct_plan->D);
  
  this_iplan->dot_r_iter_old = this_iplan->dot_r_iter;
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter =
      dot_w_complex(this_iplan->r_iter, this_iplan->w,
                    this_iplan->direct_plan->D);
  else
    this_iplan->dot_r_iter =
      dot_complex(this_iplan->r_iter, this_iplan->direct_plan->D);
  
  /** step 12
    */
  this_iplan->beta_iter =
    this_iplan->dot_r_iter / this_iplan->dot_r_iter_old;
  
  /** step 13-16
    *  updates this_iplan->f_hat_iter_2nd
    */
  if(this_iplan->infsft_flags_old & NFSFT_ITERATE_2nd)
  {
    this_iplan->gamma_iter_old = this_iplan->gamma_iter;
    this_iplan->gamma_iter =
      this_iplan->beta_iter * this_iplan->gamma_iter_old + 1;
    
    updatec_axpby_hat_old(this_iplan->f_hat_iter_2nd, 
                  this_iplan->beta_iter * this_iplan->gamma_iter_old / 
                  this_iplan->gamma_iter, this_iplan->f_hat_iter, 
                  1.0 / this_iplan->gamma_iter, this_iplan->direct_plan->M);
  }
  
  /** step 16
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->direct_plan->f_hat
    *  updates this_iplan->p_hat_iter
    */
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
    cp_w_complex(this_iplan->direct_plan->f, this_iplan->w,
            this_iplan->r_iter, this_iplan->direct_plan->D);
  else
    cp_complex(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->D); 
  
  nfsft_adjoint_old(this_iplan->direct_plan);
  
  updatec_axpy_hat_old(this_iplan->p_hat_iter, this_iplan->beta_iter,
               this_iplan->direct_plan->f_hat, this_iplan->direct_plan->M);
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
    this_iplan->dot_p_hat_iter = 
      dotproductc_w_hat_old(this_iplan->p_hat_iter, this_iplan->w_hat,
                    this_iplan->direct_plan->M);
  else
    this_iplan->dot_p_hat_iter =
      dotproductc_hat_old(this_iplan->p_hat_iter, this_iplan->direct_plan->M);
}

void infsft_loop_one_step_old(infsft_plan_old this_iplan)
{
  /*if(this_iplan->infsft_flags_old & NFSFT_LANDWEBER)
    infsft_loop_one_step_landweber(this_iplan);
  
  if(this_iplan->infsft_flags_old & NFSFT_STEEPEST_DESCENT)
    infsft_loop_one_step_steepest_descent(this_iplan);*/
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNR_E)
    infsft_loop_one_step_cgnr_e(this_iplan);
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNE_R)
    infsft_loop_one_step_cgne_r(this_iplan);
}

void infsft_finalize_old(infsft_plan_old this_iplan)
{
  int n;
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_WEIGHT)
    fftw_free(this_iplan->w);
  
  if(this_iplan->infsft_flags_old & NFSFT_PRECOMPUTE_DAMP)
  {  
    for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
    {
      free(this_iplan->w_hat[n]);
    }
    fftw_free(this_iplan->w_hat);
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_ITERATE_2nd)
  {
    for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
    {
      free(this_iplan->f_hat_iter_2nd[n]);
    }
    fftw_free(this_iplan->f_hat_iter_2nd);
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_CGNR_E)
  {
    fftw_free(this_iplan->v_iter);
    for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
    {
      free(this_iplan->z_hat_iter[n]);
    }
    fftw_free(this_iplan->z_hat_iter);
  }
  
  if(this_iplan->infsft_flags_old & NFSFT_STEEPEST_DESCENT)
    fftw_free(this_iplan->v_iter);
  
  for (n = 0; n < 2*this_iplan->direct_plan->M+1; n++)
  {
    free(this_iplan->p_hat_iter[n]);
    free(this_iplan->f_hat_iter[n]);
  }
  fftw_free(this_iplan->p_hat_iter);
  fftw_free(this_iplan->f_hat_iter);
  
  fftw_free(this_iplan->r_iter);
  fftw_free(this_iplan->given_f);
  
  free(this_iplan);
}
