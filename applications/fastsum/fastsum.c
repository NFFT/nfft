/*! \file fastsum.c
 *  \brief Fast NFFT-based summation algorithm.
 *
 *  \author Markus Fenn
 *  \date 2003-2006
 */

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "nfft3.h"
#include "fastsum.h"

/** 
 * \addtogroup applications_fastsum
 * \{
 */

/** factorial */
double fak(int n)
{
  if (n<=1) return 1.0;
  else return (double)n*fak(n-1);
}

/** binomial coefficient */
double binom(int n, int m)
{
  return fak(n)/fak(m)/fak(n-m);
}

/** basis polynomial for regularized kernel */
double BasisPoly(int m, int r, double xx)
{
  int k;
  double sum=0.0;

  for (k=0; k<=m-r; k++) {
    sum+=binom(m+k,k)*pow((xx+1.0)/2.0,(double)k);
  }
  return sum*pow((xx+1.0),(double)r)*pow(1.0-xx,(double)(m+1))/(1<<(m+1))/fak(r); /* 1<<(m+1) = 2^(m+1) */
}

/** regularized kernel with K_I arbitrary and K_B smooth to zero */
double regkern(double _Complex (*kernel)(double , int , const double *), double xx, int p, const double *param, double a, double b)
{
  int r;
  double sum=0.0;

  if (xx<-0.5)
    xx=-0.5;
  if (xx>0.5)
    xx=0.5;
  if ((xx>=-0.5+b && xx<=-a) || (xx>=a && xx<=0.5-b)) {
    return kernel(xx,0,param);
  }
  else if (xx<-0.5+b) {
    sum=(kernel(-0.5,0,param)+kernel(0.5,0,param))/2.0
        *BasisPoly(p-1,0,2.0*xx/b+(1.0-b)/b);
    for (r=0; r<p; r++) {
      sum+=pow(-b/2.0,(double)r)
          *kernel(-0.5+b,r,param)
          *BasisPoly(p-1,r,-2.0*xx/b+(b-1.0)/b);
    }
    return sum;
  }
  else if ((xx>-a) && (xx<a)) {
    for (r=0; r<p; r++) {
      sum+=pow(a,(double)r)
          *( kernel(-a,r,param)*BasisPoly(p-1,r,xx/a)
              +kernel( a,r,param)*BasisPoly(p-1,r,-xx/a)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx>0.5-b) {
    sum=(kernel(-0.5,0,param)+kernel(0.5,0,param))/2.0
        *BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b);
    for (r=0; r<p; r++) {
      sum+=pow(b/2.0,(double)r)
          *kernel(0.5-b,r,param)
          *BasisPoly(p-1,r,2.0*xx/b-(1.0-b)/b);
    }
    return sum;
  }
  return kernel(xx,0,param);
}

/** regularized kernel with K_I arbitrary and K_B periodized
 *  (used in 1D)
 */
double regkern1(double _Complex (*kernel)(double , int , const double *), double xx, int p, const double *param, double a, double b)
{
  int r;
  double sum=0.0;

  if (xx<-0.5)
    xx=-0.5;
  if (xx>0.5)
    xx=0.5;
  if ((xx>=-0.5+b && xx<=-a) || (xx>=a && xx<=0.5-b))
  {
    return kernel(xx,0,param);
  }
  else if ((xx>-a) && (xx<a))
  {
    for (r=0; r<p; r++) {
      sum+=pow(a,(double)r)
          *( kernel(-a,r,param)*BasisPoly(p-1,r,xx/a)
              +kernel( a,r,param)*BasisPoly(p-1,r,-xx/a)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx<-0.5+b)
  {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)
          *( kernel(0.5-b,r,param)*BasisPoly(p-1,r,(xx+0.5)/b)
              +kernel(-0.5+b,r,param)*BasisPoly(p-1,r,-(xx+0.5)/b)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx>0.5-b)
  {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)
          *( kernel(0.5-b,r,param)*BasisPoly(p-1,r,(xx-0.5)/b)
              +kernel(-0.5+b,r,param)*BasisPoly(p-1,r,-(xx-0.5)/b)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  return kernel(xx,0,param);
}

/** regularized kernel for even kernels with K_I even and K_B mirrored */
double regkern2(double _Complex (*kernel)(double , int , const double *), double xx, int p, const double *param, double a, double b)
{
  int r;
  double sum=0.0;

  xx=fabs(xx);

  if (xx>0.5) {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)*kernel(0.5-b,r,param)
          *(BasisPoly(p-1,r,0)+BasisPoly(p-1,r,0));
    }
    return sum;
  }
  else if ((a<=xx) && (xx<=0.5-b)) {
    return kernel(xx,0,param);
  }
  else if (xx<a) {
    for (r=0; r<p; r++) {
      sum+=pow(-a,(double)r)*kernel(a,r,param)
          *(BasisPoly(p-1,r,xx/a)+BasisPoly(p-1,r,-xx/a));
    }
    return sum;
  }
  else if ((0.5-b<xx) && (xx<=0.5)) {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)*kernel(0.5-b,r,param)
          *(BasisPoly(p-1,r,(xx-0.5)/b)+BasisPoly(p-1,r,-(xx-0.5)/b));
    }
    return sum;
  }
  return 0.0;
}

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth to K(1/2) (used in dD, d>1)
 */
double regkern3(double _Complex (*kernel)(double , int , const double *), double xx, int p, const double *param, double a, double b)
{
  int r;
  double sum=0.0;

  xx=fabs(xx);

  if (xx>=0.5) {
    /*return kern(typ,c,0,0.5);*/
    xx=0.5;
  }
  /* else */
  if ((a<=xx) && (xx<=0.5-b)) {
    return kernel(xx,0,param);
  }
  else if (xx<a) {
    for (r=0; r<p; r++) {
      sum+=pow(-a,(double)r)*kernel(a,r,param)
          *(BasisPoly(p-1,r,xx/a)+BasisPoly(p-1,r,-xx/a));
    }
    /*sum=kern(typ,c,0,xx); */
    return sum;
  }
  else if ((0.5-b<xx) && (xx<=0.5)) {
    sum=kernel(0.5,0,param)*BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b);
    /* sum=regkern2(typ,c,p,a,b, 0.5)*BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b); */
    for (r=0; r<p; r++) {
      sum+=pow(b/2.0,(double)r)
          *kernel(0.5-b,r,param)
          *BasisPoly(p-1,r,2.0*xx/b-(1.0-b)/b);
    }
    return sum;
  }
  printf("hier ");
  return 0.0;
}

/** cubic spline interpolation in near field with even kernels */
double kubintkern(double x, double *Add, int Ad, double a)
{
  double c;
  int r;
  double f0,f1,f2,f3,c1,c2,c3,c4;
  c=x*Ad/a;
  r=c; r=abs(r);
  if (r==0) {f0=Add[r+1];f1=Add[r];f2=Add[r+1];f3=Add[r+2];}
  else { f0=Add[r-1];f1=Add[r];f2=Add[r+1];f3=Add[r+2];}
  c=fabs(c);
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  c4=c1-2.0;
  /* return(-f0*(c-r)*(c-r-1.0)*(c-r-2.0)/6.0+f1*(c-r+1.0)*(c-r-1.0)*(c-r-2.0)/2-
     f2*(c-r+1.0)*(c-r)*(c-r-2.0)/2+f3*(c-r+1.0)*(c-r)*(c-r-1.0)/6.0); */
  return(-f0*c1*c3*c4/6.0+f1*c2*c3*c4/2.0-f2*c2*c1*c4/2.0+f3*c2*c1*c3/6.0);
}

/** cubic spline interpolation in near field with arbitrary kernels */
double kubintkern1(double x, double *Add, int Ad, double a)
{
  double c;
  int r;
  double f0,f1,f2,f3,c1,c2,c3,c4;
  Add+=2;
  c=(x+a)*Ad/2/a;
  r=c; r=abs(r);
  /*if (r==0) {f0=Add[r];f1=Add[r];f2=Add[r+1];f3=Add[r+2];}
  else */
  { f0=Add[r-1];f1=Add[r];f2=Add[r+1];f3=Add[r+2];}
  c=fabs(c);
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  c4=c1-2.0;
  /* return(-f0*(c-r)*(c-r-1.0)*(c-r-2.0)/6.0+f1*(c-r+1.0)*(c-r-1.0)*(c-r-2.0)/2-
     f2*(c-r+1.0)*(c-r)*(c-r-2.0)/2+f3*(c-r+1.0)*(c-r)*(c-r-1.0)/6.0); */
  return(-f0*c1*c3*c4/6.0+f1*c2*c3*c4/2.0-f2*c2*c1*c4/2.0+f3*c2*c1*c3/6.0);
}

/** quicksort algorithm for source knots and associated coefficients */
void quicksort(int d, int t, double *x, double _Complex *alpha, int N)
{
  int lpos=0;
  int rpos=N-1;
  /*double pivot=x[((N-1)/2)*d+t];*/
  double pivot=x[(N/2)*d+t];

  int k;
  double temp1;
  double _Complex temp2;

  while (lpos<=rpos)
  {
    while (x[lpos*d+t]<pivot)
      lpos++;
    while (x[rpos*d+t]>pivot)
      rpos--;
    if (lpos<=rpos)
    {
      for (k=0; k<d; k++)
      {
        temp1=x[lpos*d+k];
        x[lpos*d+k]=x[rpos*d+k];
        x[rpos*d+k]=temp1;
      }
      temp2=alpha[lpos];
      alpha[lpos]=alpha[rpos];
      alpha[rpos]=temp2;

      lpos++;
      rpos--;
    }
  }
  if (0<rpos)
    quicksort(d,t,x,alpha,rpos+1);
  if (lpos<N-1)
    quicksort(d,t,x+lpos*d,alpha+lpos,N-lpos);
}

/** recursive sort of source knots dimension by dimension to get tree structure */
void BuildTree(int d, int t, double *x, double _Complex *alpha, int N)
{
  if (N>1)
  {
    int m=N/2;

    quicksort(d,t,x,alpha,N);

    BuildTree(d, (t+1)%d, x, alpha, m);
    BuildTree(d, (t+1)%d, x+(m+1)*d, alpha+(m+1), N-m-1);
  }
}

/** fast search in tree of source knots for near field computation*/
double _Complex SearchTree(int d, int t, double *x, double _Complex *alpha, double *xmin, double *xmax, int N, double _Complex (*kernel)(double , int , const double *), const double *param, int Ad, double *Add, int p, unsigned flags)
{
  int m=N/2;
  double Min=xmin[t], Max=xmax[t], Median=x[m*d+t];
  double a=fabs(Max-Min)/2;
  int l;
  int E=0;
  double r;
  double _Complex result=0.0;

  if (N==0)
    return 0.0;

  if (Min>Median)
    result += SearchTree(d,(t+1)%d,x+(m+1)*d,alpha+(m+1),xmin,xmax,N-m-1,kernel,param,Ad,Add,p,flags);
  else if (Max<Median)
    result += SearchTree(d,(t+1)%d,x,alpha,xmin,xmax,m,kernel,param,Ad,Add,p,flags);
  else
  {
    E=0;
    for (l=0; l<d; l++)
    {
      if ( x[m*d+l]>xmin[l] && x[m*d+l]<xmax[l] )
        E++;
    }
    if (E==d)
    {
      if (d==1)
        r = xmin[0]+a-x[m];  /* remember: xmin+a = y */
      else
      {
        r=0.0;
        for (l=0; l<d; l++)
          r+=(xmin[l]+a-x[m*d+l])*(xmin[l]+a-x[m*d+l]);  /* remember: xmin+a = y */
        r=sqrt(r);
      }
      if (fabs(r)<a)
      {
        result += alpha[m]*kernel(r,0,param);                         /* alpha*(kern-regkern) */
        if (d==1)
        {
          if (flags & EXACT_NEARFIELD)
            result -= alpha[m]*regkern1(kernel,r,p,param,a,1.0/16.0); /* exact value (in 1D)  */
          else
            result -= alpha[m]*kubintkern1(r,Add,Ad,a);               /* spline approximation */
        }
        else
        {
          if (flags & EXACT_NEARFIELD)
            result -= alpha[m]*regkern(kernel,r,p,param,a,1.0/16.0);  /* exact value (in dD)  */
          else
            result -= alpha[m]*kubintkern(r,Add,Ad,a);                /* spline approximation */
        }
      }
    }
    result += SearchTree(d,(t+1)%d,x+(m+1)*d,alpha+(m+1),xmin,xmax,N-m-1,kernel,param,Ad,Add,p,flags);
    result += SearchTree(d,(t+1)%d,x,alpha,xmin,xmax,m,kernel,param,Ad,Add,p,flags);
  }
  return result;
}

/** initialization of fastsum plan */
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total, double _Complex (*kernel)(), double *param, unsigned flags, int nn, int m, int p, double eps_I, double eps_B)
{
  int t;
  int N[d], n[d];
  int n_total;

  ths->d = d;

  ths->N_total = N_total;
  ths->M_total = M_total;

  ths->x = (double *)malloc(d*N_total*(sizeof(double)));
  ths->alpha = (double _Complex *)malloc(N_total*(sizeof(double _Complex)));

  ths->y = (double *)malloc(d*M_total*(sizeof(double)));
  ths->f = (double _Complex *)malloc(M_total*(sizeof(double _Complex)));

  ths->kernel = kernel;
  ths->kernel_param = param;

  ths->flags = flags;

  ths->p = p;
  ths->eps_I = eps_I; /* =(double)ths->p/(double)nn; */  /** inner boundary */
  ths->eps_B = eps_B; /* =1.0/16.0; */                   /** outer boundary */

  /** init spline for near field computation */
  if (!(ths->flags & EXACT_NEARFIELD))
  {
    if (ths->d==1)
    {
      ths->Ad = 4*(ths->p)*(ths->p);
      ths->Add = (double *)malloc((ths->Ad+5)*(sizeof(double)));
    }
    else
    {
      ths->Ad = 2*(ths->p)*(ths->p);
      ths->Add = (double *)malloc((ths->Ad+3)*(sizeof(double)));
    }
  }

  /** init d-dimensional NFFT plan */
  ths->n = nn;
  for (t=0; t<d; t++)
  {
    N[t] = nn;
    n[t] = 2*nn;
  }
  nfft_init_guru(&(ths->mv1), d, N, N_total, n, m,
                   PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                   FFTW_MEASURE| FFTW_DESTROY_INPUT);
  nfft_init_guru(&(ths->mv2), d, N, M_total, n, m,
                   PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                   FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init d-dimensional FFTW plan */
  n_total = 1;
  for (t=0; t<d; t++)
    n_total *= nn;

  ths->b = (fftw_complex *)fftw_malloc(n_total*sizeof(fftw_complex));
  ths->fft_plan = fftw_plan_dft(d,N,ths->b,ths->b,FFTW_FORWARD,FFTW_ESTIMATE);

}

/** finalization of fastsum plan */
void fastsum_finalize(fastsum_plan *ths)
{
  free(ths->x);
  free(ths->alpha);
  free(ths->y);
  free(ths->f);

  if (!(ths->flags & EXACT_NEARFIELD))
    free(ths->Add);

  nfft_finalize(&(ths->mv1));
  nfft_finalize(&(ths->mv2));

  fftw_destroy_plan(ths->fft_plan);
  fftw_free(ths->b);
}

/** direct computation of sums */
void fastsum_exact(fastsum_plan *ths)
{
  int j,k;
  int t;
  double r;

  for (j=0; j<ths->M_total; j++)
  {
    ths->f[j]=0.0;
    for (k=0; k<ths->N_total; k++)
    {
      if (ths->d==1)
        r = ths->y[j] - ths->x[k];
      else
      {
        r=0.0;
        for (t=0; t<ths->d; t++)
          r += (ths->y[j*ths->d+t]-ths->x[k*ths->d+t])*(ths->y[j*ths->d+t]-ths->x[k*ths->d+t]);
        r=sqrt(r);
      }
      ths->f[j] += ths->alpha[k] * ths->kernel(r,0,ths->kernel_param);
    }
  }
}

/** precomputation for fastsum */
void fastsum_precompute(fastsum_plan *ths)
{
  int j,k,t;
  int n_total;

  /** sort source knots */
  BuildTree(ths->d,0,ths->x,ths->alpha,ths->N_total);

  /** precompute spline values for near field*/
  if (!(ths->flags & EXACT_NEARFIELD))
  {
    if (ths->d==1)
      for (k=-ths->Ad/2-2; k <= ths->Ad/2+2; k++)
        ths->Add[k+ths->Ad/2+2] = regkern1(ths->kernel, ths->eps_I*(double)k/ths->Ad*2, ths->p, ths->kernel_param, ths->eps_I, ths->eps_B);
    else
      for (k=0; k <= ths->Ad+2; k++)
        ths->Add[k] = regkern3(ths->kernel, ths->eps_I*(double)k/ths->Ad, ths->p, ths->kernel_param, ths->eps_I, ths->eps_B);
  }

  /** init NFFT plan for transposed transform in first step*/
  for (k=0; k<ths->mv1.M_total; k++)
    for (t=0; t<ths->mv1.d; t++)
      ths->mv1.x[ths->mv1.d*k+t] = - ths->x[ths->mv1.d*k+t];  /* note the factor -1 for transposed transform instead of adjoint*/

  /** precompute psi, the entries of the matrix B */
  if(ths->mv1.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&(ths->mv1));

  if(ths->mv1.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&(ths->mv1));

  if(ths->mv1.nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&(ths->mv1));

  /** init Fourier coefficients */
  for(k=0; k<ths->mv1.M_total;k++)
    ths->mv1.f[k] = ths->alpha[k];

  /** init NFFT plan for transform in third step*/
  for (j=0; j<ths->mv2.M_total; j++)
    for (t=0; t<ths->mv2.d; t++)
      ths->mv2.x[ths->mv2.d*j+t] = - ths->y[ths->mv2.d*j+t];  /* note the factor -1 for conjugated transform instead of standard*/

  /** precompute psi, the entries of the matrix B */
  if(ths->mv2.nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&(ths->mv2));

  if(ths->mv2.nfft_flags & PRE_PSI)
    nfft_precompute_psi(&(ths->mv2));

  if(ths->mv2.nfft_flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&(ths->mv2));


  /** precompute Fourier coefficients of regularised kernel*/
  n_total = 1;
  for (t=0; t<ths->d; t++)
    n_total *= ths->n;

  for (j=0; j<n_total; j++)
  {
    if (ths->d==1)
      ths->b[j] = regkern1(ths->kernel, (double)j / (ths->n) - 0.5, ths->p, ths->kernel_param, ths->eps_I, ths->eps_B)/n_total;
    else
    {
      k=j;
      ths->b[j]=0.0;
      for (t=0; t<ths->d; t++)
      {
        ths->b[j] += ((double)(k % (ths->n)) / (ths->n) - 0.5) * ((double)(k % (ths->n)) / (ths->n) - 0.5);
        k = k / (ths->n);
      }
      ths->b[j] = regkern3(ths->kernel, sqrt(ths->b[j]), ths->p, ths->kernel_param, ths->eps_I, ths->eps_B)/n_total;
    }
  }

  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);
  fftw_execute(ths->fft_plan);
  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);

}

/** fast NFFT-based summation */
void fastsum_trafo(fastsum_plan *ths)
{
  int j,k,t;
  double *ymin, *ymax;   /** limits for d-dimensional near field box */

  ymin = (double *)malloc(ths->d*(sizeof(double)));
  ymax = (double *)malloc(ths->d*(sizeof(double)));

  /** first step of algorithm */
  nfft_adjoint(&(ths->mv1));

  /** second step of algorithm */
  for (k=0; k<ths->mv2.N_total; k++)
    ths->mv2.f_hat[k] = ths->b[k] * ths->mv1.f_hat[k];

  /** third step of algorithm */
  nfft_trafo(&(ths->mv2));

  /** add near field */
  for (j=0; j<ths->M_total; j++)
  {
    for (t=0; t<ths->d; t++)
    {
      ymin[t] = ths->y[ths->d*j+t] - ths->eps_I;
      ymax[t] = ths->y[ths->d*j+t] + ths->eps_I;
    }
    ths->f[j] = ths->mv2.f[j] + SearchTree(ths->d,0, ths->x, ths->alpha, ymin, ymax, ths->N_total, ths->kernel, ths->kernel_param, ths->Ad, ths->Add, ths->p, ths->flags);
    /* ths->f[j] = ths->mv2.f[j]; */
    /* ths->f[j] = SearchTree(ths->d,0, ths->x, ths->alpha, ymin, ymax, ths->N_total, ths->kernel, ths->kernel_param, ths->Ad, ths->Add, ths->p, ths->flags); */
  }
}
/* \} */

/* fastsum.c */
