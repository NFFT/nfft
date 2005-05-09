/** inverse nfft 
* -----------------------------------------------------------------------------
* -----------------------------------------------------------------------------
*/
/*void infft_init_specific(infft_plan *this_iplan, nfft_plan *direct_plan,
                         int infft_flags)
{
  this_iplan->direct_plan = direct_plan;
  this_iplan->infft_flags = infft_flags;
  
  this_iplan->given_f = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));
  
  this_iplan->r_iter = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));
  
  this_iplan->f_hat_iter = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));
  
  this_iplan->p_hat_iter = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));
  
  if(this_iplan->infft_flags & LANDWEBER)
  {
    this_iplan->z_hat_iter = this_iplan->p_hat_iter;
  }
  
  if(this_iplan->infft_flags & STEEPEST_DESCENT)
  {
    this_iplan->z_hat_iter = this_iplan->p_hat_iter;
    
    this_iplan->v_iter = (fftw_complex*)
      fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));
  }
  
  if(this_iplan->infft_flags & CGNR_E)
  {
    this_iplan->z_hat_iter = (fftw_complex*)
    fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));
    
    this_iplan->v_iter = (fftw_complex*)
      fftw_malloc(this_iplan->direct_plan->M*sizeof(fftw_complex));
  }
  
  if(this_iplan->infft_flags & CGNE_R)
  {
    this_iplan->z_hat_iter = this_iplan->p_hat_iter;
  }
  
  if(this_iplan->infft_flags & ITERATE_2nd)
    this_iplan->f_hat_iter_2nd = (fftw_complex*) 
      fftw_malloc(this_iplan->direct_plan->N_L*sizeof(fftw_complex));
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->w = 
      (double*) fftw_malloc(this_iplan->direct_plan->M*sizeof(double));
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->w_hat = 
      (double*) fftw_malloc(this_iplan->direct_plan->N_L*sizeof(double));
}*/

/*void infft_init(infft_plan *this_iplan, nfft_plan *direct_plan)
{
  infft_init_specific(this_iplan, direct_plan, CGNR_E);
}*/

/*void infft_before_loop_help(infft_plan *this_iplan)
{*/
  /** step 2
  *  overwrites this_iplan->direct_plan->f_hat 
  *  overwrites this_iplan->r_iter
  */
  /*copyc(this_iplan->direct_plan->f_hat, this_iplan->f_hat_iter,
        this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->r_iter, this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->r_iter, this_iplan->direct_plan->f);
  
  updatec_axpy(this_iplan->r_iter, -1.0, this_iplan->given_f,
               this_iplan->direct_plan->M);
  
  if((!(this_iplan->infft_flags & LANDWEBER)) ||
     (this_iplan->infft_flags & NORMS_FOR_LANDWEBER))
  {
    if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
      this_iplan->dot_r_iter =
        dotproductc_w(this_iplan->r_iter, this_iplan->w,
                      this_iplan->direct_plan->M);
    else
      this_iplan->dot_r_iter =
        dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);
  }*/
  
  /** step 3
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->z_hat_iter resp. this_iplan->z_hat_iter
    */ 
  
  /*if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w, this_iplan->r_iter,
            this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter, this_iplan->direct_plan->f_hat);
  
  if((!(this_iplan->infft_flags & LANDWEBER)) ||
     (this_iplan->infft_flags & NORMS_FOR_LANDWEBER))
  {
    if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
      this_iplan->dot_z_hat_iter = 
        dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat,
                      this_iplan->direct_plan->N_L);
    else
      this_iplan->dot_z_hat_iter =
        dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  }
  
  if(this_iplan->infft_flags & CGNE_R)
    this_iplan->dot_p_hat_iter = this_iplan->dot_z_hat_iter;
  
}*/ /* void infft_before_loop_help */

/*void infft_before_loop(infft_plan *this_iplan)
{
  infft_before_loop_help(this_iplan);
  
  if(this_iplan->infft_flags & CGNR_E)
  {*/
    /** step 4-6
    *  overwrites this_iplan->f_hat_iter_2nd
    */
/*    if(this_iplan->infft_flags & ITERATE_2nd)
      copyc(this_iplan->f_hat_iter_2nd, this_iplan->f_hat_iter,
            this_iplan->direct_plan->N_L);
    
    /** step 7
    *  overwrites this_iplan->p_hat_iter
    */
    /*copyc(this_iplan->p_hat_iter, this_iplan->z_hat_iter,
          this_iplan->direct_plan->N_L);
  }
  
  if(this_iplan->infft_flags & CGNE_R)
  {*/
    /** step 4-7
    *  overwrites this_iplan->f_hat_iter_2nd
    */
    /*if(this_iplan->infft_flags & ITERATE_2nd)
    {
      this_iplan->gamma_iter=1.0;
      copyc(this_iplan->f_hat_iter_2nd, this_iplan->f_hat_iter,
            this_iplan->direct_plan->N_L);
    }
  }
}*/ /* void infft_before_loop */

/*void infft_loop_one_step_landweber(infft_plan *this_iplan)
{*/
  /** step 5
  *  updates this_iplan->f_hat_iter
  */
  /*if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
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
  
  if(this_iplan->infft_flags & NORMS_FOR_LANDWEBER)
  {
    if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
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
  /*if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w,
            this_iplan->r_iter, this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  
  if(this_iplan->infft_flags & NORMS_FOR_LANDWEBER)
  {
    if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
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
  /*if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    copyc_w(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
            this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  else
    copyc(this_iplan->direct_plan->f_hat, this_iplan->z_hat_iter,
          this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
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
  /*if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
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
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter = dotproductc_w(this_iplan->r_iter, this_iplan->w,
                                           this_iplan->direct_plan->M);
  else
    this_iplan->dot_r_iter =
      dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);*/
  
  /** step 9
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->z_hat_iter
    */
  /*if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w, this_iplan->r_iter,
            this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->dot_z_hat_iter =
      dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat, 
                    this_iplan->direct_plan->N_L);
  else
    this_iplan->dot_z_hat_iter =
      dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
}*/ /* void infft_loop_one_step_steepest_descent */

/*void infft_loop_one_step_cgnr_e(infft_plan *this_iplan)
{*/
  /** step 9
  *  overwrites this_iplan->direct_plan->f_hat 
  *  overwrites this_iplan->v_iter
  */
  /*if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    copyc_w(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
            this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);
  else
    copyc(this_iplan->direct_plan->f_hat, this_iplan->p_hat_iter,
          this_iplan->direct_plan->N_L);
  
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  nfft_trafo(this_iplan->direct_plan);
  SWAPC(this_iplan->v_iter,this_iplan->direct_plan->f);
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_v_iter = dotproductc_w(this_iplan->v_iter, this_iplan->w,
                                           this_iplan->direct_plan->M);
  else
    this_iplan->dot_v_iter =
      dotproductc(this_iplan->v_iter, this_iplan->direct_plan->M);*/
  
  /** step 10
    */
  //this_iplan->alpha_iter =
  //  this_iplan->dot_z_hat_iter / this_iplan->dot_v_iter;
  
  /** step 11
    *  updates this_iplan->f_hat_iter
    */
  /*if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                  this_iplan->w_hat, this_iplan->p_hat_iter,
                  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter, 
                 this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);*/
  
  /** step 12-15
    *  updates this_iplan->f_hat_iter_2nd
    */
  /*if(this_iplan->infft_flags & ITERATE_2nd)
  {
    this_iplan->alpha_iter_2nd =
    this_iplan->dot_r_iter / this_iplan->dot_z_hat_iter;
    
    if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
      updatec_xpawy(this_iplan->f_hat_iter_2nd, this_iplan->alpha_iter_2nd,
                    this_iplan->w_hat, this_iplan->z_hat_iter,
                    this_iplan->direct_plan->N_L);
    else
      updatec_xpay(this_iplan->f_hat_iter_2nd, this_iplan->alpha_iter_2nd,
                   this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
  }*/
  
  /** step 16
    *  updates this_iplan->r_iter
    */
  /*updatec_xpay(this_iplan->r_iter, -this_iplan->alpha_iter, this_iplan->v_iter,
               this_iplan->direct_plan->M);
  
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter = dotproductc_w(this_iplan->r_iter, this_iplan->w,
                                           this_iplan->direct_plan->M);
  else
    this_iplan->dot_r_iter =
      dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);*/
  
  /** step 17
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->direct_plan->r_iter
    */
  /*if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w, 
            this_iplan->r_iter, this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->M);
  
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  nfft_adjoint(this_iplan->direct_plan);
  SWAPC(this_iplan->z_hat_iter,this_iplan->direct_plan->f_hat);
  
  this_iplan->dot_z_hat_iter_old = this_iplan->dot_z_hat_iter;
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->dot_z_hat_iter =
      dotproductc_w(this_iplan->z_hat_iter, this_iplan->w_hat, 
                    this_iplan->direct_plan->N_L);
  else
    this_iplan->dot_z_hat_iter =
      dotproductc(this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);*/
  
  /** step 18
    */
  //this_iplan->beta_iter =
  //  this_iplan->dot_z_hat_iter / this_iplan->dot_z_hat_iter_old;
  
  /** step 19
    *  updates this_iplan->p_hat_iter
    */
  /*updatec_axpy(this_iplan->p_hat_iter, this_iplan->beta_iter, 
               this_iplan->z_hat_iter, this_iplan->direct_plan->N_L);
}*/ /* void infft_loop_one_step_cgnr_e */

/*void infft_loop_one_step_cgne_r(infft_plan *this_iplan)
{*/
  /** step 9
  */
//  this_iplan->alpha_iter =
//  this_iplan->dot_r_iter / this_iplan->dot_p_hat_iter;
  
  /** step 10
  *  updates this_iplan->f_hat_iter
  */
/*  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    updatec_xpawy(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                  this_iplan->w_hat, this_iplan->p_hat_iter,
                  this_iplan->direct_plan->N_L);
  else
    updatec_xpay(this_iplan->f_hat_iter, this_iplan->alpha_iter,
                 this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);*/
  
  /** step 11
    *  overwrites this_iplan->direct_plan->f_hat 
    *  overwrites this_iplan->direct_plan->f
    *  updates this_iplan->r_iter
    */
/*  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    copyc_w(this_iplan->direct_plan->f_hat, this_iplan->w_hat,
            this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);
  else
    copyc(this_iplan->direct_plan->f_hat, this_iplan->p_hat_iter,
          this_iplan->direct_plan->N_L);
  
  nfft_trafo(this_iplan->direct_plan);
  
  updatec_xpay(this_iplan->r_iter, -this_iplan->alpha_iter,
               this_iplan->direct_plan->f, this_iplan->direct_plan->M);
  
  this_iplan->dot_r_iter_old = this_iplan->dot_r_iter;
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    this_iplan->dot_r_iter =
      dotproductc_w(this_iplan->r_iter, this_iplan->w,
                    this_iplan->direct_plan->M);
  else
    this_iplan->dot_r_iter =
      dotproductc(this_iplan->r_iter, this_iplan->direct_plan->M);*/
  
  /** step 12
    */
//  this_iplan->beta_iter =
//    this_iplan->dot_r_iter / this_iplan->dot_r_iter_old;
  
  /** step 13-16
    *  updates this_iplan->f_hat_iter_2nd
    */
/*  if(this_iplan->infft_flags & ITERATE_2nd)
  {
    this_iplan->gamma_iter_old = this_iplan->gamma_iter;
    this_iplan->gamma_iter =
      this_iplan->beta_iter * this_iplan->gamma_iter_old + 1;
    
    updatec_axpby(this_iplan->f_hat_iter_2nd, 
                  this_iplan->beta_iter * this_iplan->gamma_iter_old / 
                  this_iplan->gamma_iter, this_iplan->f_hat_iter, 
                  1.0 / this_iplan->gamma_iter, this_iplan->direct_plan->N_L);
  }*/
  
  /** step 16
    *  overwrites this_iplan->direct_plan->f
    *  overwrites this_iplan->direct_plan->f_hat
    *  updates this_iplan->p_hat_iter
    */
  /*if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    copyc_w(this_iplan->direct_plan->f, this_iplan->w,
            this_iplan->r_iter, this_iplan->direct_plan->M);
  else
    copyc(this_iplan->direct_plan->f, this_iplan->r_iter,
          this_iplan->direct_plan->M); 
  
  nfft_adjoint(this_iplan->direct_plan);
  
  updatec_axpy(this_iplan->p_hat_iter, this_iplan->beta_iter,
               this_iplan->direct_plan->f_hat, this_iplan->direct_plan->N_L);
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    this_iplan->dot_p_hat_iter = 
      dotproductc_w(this_iplan->p_hat_iter, this_iplan->w_hat,
                    this_iplan->direct_plan->N_L);
  else
    this_iplan->dot_p_hat_iter =
      dotproductc(this_iplan->p_hat_iter, this_iplan->direct_plan->N_L);
}

void infft_loop_one_step(infft_plan *this_iplan)
{
  if(this_iplan->infft_flags & LANDWEBER)
    infft_loop_one_step_landweber(this_iplan);
  
  if(this_iplan->infft_flags & STEEPEST_DESCENT)
    infft_loop_one_step_steepest_descent(this_iplan);
  
  if(this_iplan->infft_flags & CGNR_E)
    infft_loop_one_step_cgnr_e(this_iplan);
  
  if(this_iplan->infft_flags & CGNE_R)
    infft_loop_one_step_cgne_r(this_iplan);
}

void infft_finalize(infft_plan *this_iplan)
{
  if(this_iplan->infft_flags & PRECOMPUTE_WEIGHT)
    fftw_free(this_iplan->w);
  
  if(this_iplan->infft_flags & PRECOMPUTE_DAMP)
    fftw_free(this_iplan->w_hat);
  
  if(this_iplan->infft_flags & ITERATE_2nd)
    fftw_free(this_iplan->f_hat_iter_2nd);
  
  if(this_iplan->infft_flags & CGNR_E)
  {
    fftw_free(this_iplan->v_iter);
    fftw_free(this_iplan->z_hat_iter);
  }
  
  if(this_iplan->infft_flags & STEEPEST_DESCENT)
    fftw_free(this_iplan->v_iter);
  
  fftw_free(this_iplan->p_hat_iter);
  fftw_free(this_iplan->f_hat_iter);
  
  fftw_free(this_iplan->r_iter);
  fftw_free(this_iplan->given_f);
}*/
