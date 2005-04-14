  
  int r, k1, k2, a, b;
  const int N=pow(2.0,j);

  /* clear memory */
  memset(f,0,N*N*sizeof(fftw_complex));

  /* set random values at hyperbolic grid point */
  for (r=0; r<ceil(j/2.0); r++){
    a=pow(2,j-r-2); b=pow(2,r);
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        /* top */
        f[( k1+a+N/2 )*N+(k2-b/2+N/2)][0]=(double)rand()/RAND_MAX;
        f[( k1+a+N/2 )*N+(k2-b/2+N/2)][1]=(double)rand()/RAND_MAX;
        /* bottom */
        f[(k1-2*a+N/2)*N+(k2-b/2+N/2)][0]=(double)rand()/RAND_MAX;
        f[(k1-2*a+N/2)*N+(k2-b/2+N/2)][1]=(double)rand()/RAND_MAX;
        /* right */
        f[(k2-b/2+N/2)*N+( k1+a+N/2 )][0]=(double)rand()/RAND_MAX;
        f[(k2-b/2+N/2)*N+( k1+a+N/2 )][1]=(double)rand()/RAND_MAX;
        /* left */
        f[(k2-b/2+N/2)*N+(k1-2*a+N/2)][0]=(double)rand()/RAND_MAX;
        f[(k2-b/2+N/2)*N+(k1-2*a+N/2)][1]=(double)rand()/RAND_MAX;
      }
    }
  }
  /* center */
  a=pow(2,floor(j/2));
  for (k1=0; k1<a; k1++){
    for (k2=0; k2<a; k2++){
      f[(k1-a/2+N/2)*N+(k2-a/2+N/2)][0]=(double)rand()/RAND_MAX;
      f[(k1-a/2+N/2)*N+(k2-a/2+N/2)][1]=(double)rand()/RAND_MAX;
    }
  }







  

  /* initialise nfft_plans */
  for (r=0; r<ceil(j/2.0); r++) {
    a=pow(2,j-r-2); b=pow(2,r);
    /* top */
    N[0]=a; n[0]=2*N[0];
    N[1]=b; n[1]=2*N[1];
    nfft_init_specific(&(this_sparse_plan->this_nfft_plans[4*r+0]),2,N,M,n, m, PRE_PHI_HUT | PRE_PSI, FFTW_ESTIMATE);
    /* bottom */
    N[0]=a; n[0]=2*N[0];
    N[1]=b; n[1]=2*N[1];
    nfft_init_specific(&(this_sparse_plan->this_nfft_plans[4*r+1]),2,N,M,n, m, PRE_PHI_HUT | PRE_PSI, FFTW_ESTIMATE);
    /* right */
    N[0]=b; n[0]=2*N[0];
    N[1]=a; n[1]=2*N[1];
    nfft_init_specific(&(this_sparse_plan->this_nfft_plans[4*r+2]),2,N,M,n, m, PRE_PHI_HUT | PRE_PSI, FFTW_ESTIMATE);

    /* left */
    N[0]=b; n[0]=2*N[0];
    N[1]=a; n[1]=2*N[1];
    nfft_init_specific(&(this_sparse_plan->this_nfft_plans[4*r+3]),2,N,M,n, m, PRE_PHI_HUT | PRE_PSI, FFTW_ESTIMATE);
  }

  /* center */
  a=pow(2,floor(j/2));
  N[0]=a; n[0]=2*N[0];
  N[1]=a; n[1]=2*N[1];
  nfft_init_specific(&(this_sparse_plan->this_nfft_plans[4*(int)ceil(j/2.0)]),2,N,M,n, m, PRE_PHI_HUT | PRE_PSI, FFTW_ESTIMATE);





/*========================================================*/
void sparse_nfft_copy_plan(sparse_nfft_plan *my_sparse_plan, nfft_plan *my_nfft_plan){
  int r, k1, k2, a, b;
  const int j=my_sparse_plan->j;
  const int N=pow(2.0,j);

  #define F_HAT my_nfft_plan->f_hat
  #define spNFFT my_sparse_plan->my_nfft_plans

  /* copy coefficients f_hat */
  for (r=0; r<ceil(j/2.0); r++){
    a=pow(2,j-r-2); b=pow(2,r);
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        /* top */
        spNFFT[4*r+0].f_hat[k1*b+k2][0]=F_HAT[( k1+a+N/2 )*N+(k2-b/2+N/2)][0];
        spNFFT[4*r+0].f_hat[k1*b+k2][1]=F_HAT[( k1+a+N/2 )*N+(k2-b/2+N/2)][1];
        /* bottom */
        spNFFT[4*r+1].f_hat[k1*b+k2][0]=F_HAT[(k1-2*a+N/2)*N+(k2-b/2+N/2)][0];
        spNFFT[4*r+1].f_hat[k1*b+k2][1]=F_HAT[(k1-2*a+N/2)*N+(k2-b/2+N/2)][1];
        /* right */
        spNFFT[4*r+2].f_hat[k2*a+k1][0]=F_HAT[(k2-b/2+N/2)*N+( k1+a+N/2 )][0];
        spNFFT[4*r+2].f_hat[k2*a+k1][1]=F_HAT[(k2-b/2+N/2)*N+( k1+a+N/2 )][1];
        /* left */
        spNFFT[4*r+3].f_hat[k2*a+k1][0]=F_HAT[(k2-b/2+N/2)*N+(k1-2*a+N/2)][0];
        spNFFT[4*r+3].f_hat[k2*a+k1][1]=F_HAT[(k2-b/2+N/2)*N+(k1-2*a+N/2)][1];
      }
    }
  }
  /* center */
  a=pow(2,floor(j/2));
  for (k1=0; k1<a; k1++){
    for (k2=0; k2<a; k2++){
      spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][0]=F_HAT[(k1-a/2+N/2)*N+(k2-a/2+N/2)][0];
      spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][1]=F_HAT[(k1-a/2+N/2)*N+(k2-a/2+N/2)][1];
    }
  }

  /* copy knots x */
  for (r=0; r<=4*ceil(j/2.0); r++) {
    memcpy(spNFFT[r].x, my_nfft_plan->x, (my_sparse_plan->M)*sizeof(fftw_complex));
    //free(spNFFT[r].x);
    //spNFFT[r].x=my_nfft_plan->x;
  }

  /* precompute psi */
  if(my_nfft_plan->nfft_flags & PRE_PSI) {
    for (r=0; r<=4*ceil(j/2.0); r++) {
      nfft_precompute_psi(&(spNFFT[r]));
    }
  }

  /* premultiply factors to psi */

  #undef F_HAT
  #undef spNFFT
}



/*========================================================*/
void sparse_nfft_compare(sparse_nfft_plan *my_sparse_plan, nfft_plan *my_nfft_plan)
{
  int r, k1, k2, a, b;
  const int j=my_sparse_plan->j;
  const int N=pow(2.0,j);
  const int M=my_sparse_plan->M;

  #define F_HAT my_nfft_plan->f_hat
  #define spNFFT my_sparse_plan->my_nfft_plans

  /* compare coefficients f_hat */
  for (r=0; r<ceil(j/2.0); r++){
    a=pow(2,j-r-2); b=pow(2,r);
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        /* top */
        if (spNFFT[4*r+0].f_hat[k1*b+k2][0]!=F_HAT[( k1+a+N/2 )*N+(k2-b/2+N/2)][0]) printf("!");
        if (spNFFT[4*r+0].f_hat[k1*b+k2][1]!=F_HAT[( k1+a+N/2 )*N+(k2-b/2+N/2)][1]) printf("!");
        /* bottom */
        if (spNFFT[4*r+1].f_hat[k1*b+k2][0]!=F_HAT[(k1-2*a+N/2)*N+(k2-b/2+N/2)][0]) printf("!");
        if (spNFFT[4*r+1].f_hat[k1*b+k2][1]!=F_HAT[(k1-2*a+N/2)*N+(k2-b/2+N/2)][1]) printf("!");
        /* right */
        if (spNFFT[4*r+2].f_hat[k2*a+k1][0]!=F_HAT[(k2-b/2+N/2)*N+( k1+a+N/2 )][0]) printf("!");
        if (spNFFT[4*r+2].f_hat[k2*a+k1][1]!=F_HAT[(k2-b/2+N/2)*N+( k1+a+N/2 )][1]) printf("!");
        /* left */
        if (spNFFT[4*r+3].f_hat[k2*a+k1][0]!=F_HAT[(k2-b/2+N/2)*N+(k1-2*a+N/2)][0]) printf("!");
        if (spNFFT[4*r+3].f_hat[k2*a+k1][1]!=F_HAT[(k2-b/2+N/2)*N+(k1-2*a+N/2)][1]) printf("!");
      }
    }
  }
  /* center */
  a=pow(2,floor(j/2));
  for (k1=0; k1<a; k1++){
    for (k2=0; k2<a; k2++){
      if (spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][0]!=F_HAT[(k1-a/2+N/2)*N+(k2-a/2+N/2)][0]) printf("!");
      if (spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][1]!=F_HAT[(k1-a/2+N/2)*N+(k2-a/2+N/2)][1]) printf("!");
    }
  }

  /* compare knots x */
  for (r=0; r<=4*ceil(j/2.0); r++) {
    for (k1=0; k1<M; k1++){
      if (spNFFT[r].x[2*k1+0]!=my_nfft_plan->x[2*k1+0]) printf("!");
      if (spNFFT[r].x[2*k1+1]!=my_nfft_plan->x[2*k1+1]) printf("!");
    }
  }

  #undef F_HAT
  #undef spNFFT
}




/*========================================================*/
void sparse_ndft(sparse_nfft_plan *my_sparse_plan){
  int l, r, k1, k2, a, b;
  const int j=my_sparse_plan->j;
  const int M=my_sparse_plan->M;
  double temp;

  #define spF my_sparse_plan->f
  #define spNFFT my_sparse_plan->my_nfft_plans

  //printf("M=%d\n",M);

  for (l=0; l<M; l++){
    spF[l][0]=0.0;
    spF[l][1]=0.0;
    double x=spNFFT[0].x[2*l+0];
    double y=spNFFT[0].x[2*l+1];
    //if (x!=my_plan->x[2*l+0]) printf("x");
    //if (y!=my_plan->x[2*l+1]) printf("y");
    //printf("(%e, %e)\n",x,y);
    for (r=0; r<ceil(j/2.0); r++){
      a=pow(2,j-r-2); b=pow(2,r);
      for (k1=0; k1<a; k1++){
        for (k2=0; k2<b; k2++){
          //if (spNFFT[4*r+2].f_hat[k2*a+k1][1] != my_plan->f_hat[(k2-(int)floor(b/2)+N/2)*N+( k1+a+N/2 )][1])
          //  printf("!");
          /* top */
          temp=2.0*PI*( ( k1+a )*x + (k2-(int)floor(b/2))*y );
          spF[l][0]+=spNFFT[4*r+0].f_hat[k1*b+k2][0]*cos(temp)+spNFFT[4*r+0].f_hat[k1*b+k2][1]*sin(temp);
          spF[l][1]+=spNFFT[4*r+0].f_hat[k1*b+k2][1]*cos(temp)-spNFFT[4*r+0].f_hat[k1*b+k2][0]*sin(temp);
          /* bottom */
          temp=2.0*PI*( (k1-2*a)*x + (k2-(int)floor(b/2))*y );
          spF[l][0]+=spNFFT[4*r+1].f_hat[k1*b+k2][0]*cos(temp)+spNFFT[4*r+1].f_hat[k1*b+k2][1]*sin(temp);
          spF[l][1]+=spNFFT[4*r+1].f_hat[k1*b+k2][1]*cos(temp)-spNFFT[4*r+1].f_hat[k1*b+k2][0]*sin(temp);
          /* right */
          temp=2.0*PI*( (k2-(int)floor(b/2))*x + ( k1+a )*y );
          spF[l][0]+=spNFFT[4*r+2].f_hat[k2*a+k1][0]*cos(temp)+spNFFT[4*r+2].f_hat[k2*a+k1][1]*sin(temp);
          spF[l][1]+=spNFFT[4*r+2].f_hat[k2*a+k1][1]*cos(temp)-spNFFT[4*r+2].f_hat[k2*a+k1][0]*sin(temp);
          /* left */
          temp=2.0*PI*( (k2-(int)floor(b/2))*x + (k1-2*a)*y );
          spF[l][0]+=spNFFT[4*r+3].f_hat[k2*a+k1][0]*cos(temp)+spNFFT[4*r+3].f_hat[k2*a+k1][1]*sin(temp);
          spF[l][1]+=spNFFT[4*r+3].f_hat[k2*a+k1][1]*cos(temp)-spNFFT[4*r+3].f_hat[k2*a+k1][0]*sin(temp);
        }
      }
    }
    /* center */
    a=pow(2,floor(j/2));
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<a; k2++){
        //if (spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][1] != my_plan->f_hat[(k1-(int)floor(a/2)+N/2)*N+(k2-(int)floor(a/2)+N/2)][1])
        //  printf("!");
        temp=2.0*PI*( (k1-(int)floor(a/2))*x + (k2-(int)floor(a/2))*y );
        spF[l][0]+=spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][0]*cos(temp)+spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][1]*sin(temp);
        spF[l][1]+=spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][1]*cos(temp)-spNFFT[4*(int)ceil(j/2.0)].f_hat[k1*a+k2][0]*sin(temp);
      }
    }
  }

  #undef spF
  #undef spNFFT
}

/*========================================================*/
void print_sparse_koeffs(fftw_complex* f, int N){
  int k1, k2;

  //vpr_c(f,N,"ndft");
  for (k1=0; k1<N; k1++) {
    for (k2=0; k2<N; k2++) {
      if (f[k1*N+k2][0]!=0.0)
        printf("*");
      else printf(" ");
    }
    printf("\n");
  }
}




  /* full NFFT */
  printf("(fast full   ) elapsed time: "); fflush(stdout);
  t=second();
    nfft_trafo(&my_plan);
  t=second()-t;
  printf("%.3f secs.\n",t);

  /* sparse NDFT */
  printf("(slow sparse ) elapsed time: "); fflush(stdout);
  t=second();
  //sparse_ndft(&my_sparse_plan);
    //print_sparse_koeffs(my_plan.f_hat,my_plan.N[0]);
  t=second()-t;
  printf("%.3f secs.\n",t);

  /* copy result */
  for (k=0; k<M_sparse; k++){
    slow_sparse[k][0]=my_sparse_plan.f[k][0]; my_sparse_plan.f[k][0]=0.0;
    slow_sparse[k][1]=my_sparse_plan.f[k][1]; my_sparse_plan.f[k][1]=0.0;
  }

  /* sparse NFFT */
  printf("(fast sparse ) elapsed time: "); fflush(stdout);
  t=second();
  //sparse_nfft(&my_sparse_plan);
    //print_sparse_koeffs(my_plan.f_hat,my_plan.N[0]);
  t=second()-t;
  printf("%.3f secs.\n",t);

  /* results */
  //printf("Error (slow full  -fast full  ): %e\n",fehler(slow_full  , my_plan.f       , my_plan.M));
  //printf("Error (slow sparse-fast full  ): %e\n",fehler(slow_sparse, my_plan.f       , M_sparse));
  printf("Error (slow sparse-fast full  ): %e\n",E_2_error_c(slow_sparse, my_plan.f       , M_sparse));
  //printf("Error (slow full  -slow sparse): %e\n",fehler(slow_full  , slow_sparse     , my_plan.M));
  //printf("Error (slow full  -fast sparse): %e\n",fehler(slow_full  , my_sparse_plan.f, my_plan.M));
  //printf("Error (slow sparse-fast sparse): %e\n",fehler(slow_sparse, my_sparse_plan.f, M_sparse));
  printf("Error (slow sparse-fast sparse): %e\n",E_2_error_c(slow_sparse, my_sparse_plan.f, M_sparse));
  //printf("Error: %e\n",E_2_error_c(my_plan.f,my_sparse_plan.f,my_plan.M));
