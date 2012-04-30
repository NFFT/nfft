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

/*! \file fastsum.c
 *  \brief Fast NFFT-based summation algorithm.
 *
 *  \author Markus Fenn
 *  \date 2003-2006
 */
#include "config.h"

#include <stdlib.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3util.h"
#include "nfft3.h"
#include "fastsum.h"
#include "infft.h"

/** Required for test if (ths->k == one_over_x) */
#include "kernels.h"

/**
 * \addtogroup applications_fastsum
 * \{
 */

/** max */
int max_i(int a, int b)
{
  return a >= b ? a : b;
}

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
double _Complex regkern(kernel k, double xx, int p, const double *param, double a, double b)
{
  int r;
  double _Complex sum=0.0;

  if (xx<-0.5)
    xx=-0.5;
  if (xx>0.5)
    xx=0.5;
  if ((xx>=-0.5+b && xx<=-a) || (xx>=a && xx<=0.5-b)) {
    return k(xx,0,param);
  }
  else if (xx<-0.5+b) {
    sum=(k(-0.5,0,param)+k(0.5,0,param))/2.0
        *BasisPoly(p-1,0,2.0*xx/b+(1.0-b)/b);
    for (r=0; r<p; r++) {
      sum+=pow(-b/2.0,(double)r)
          *k(-0.5+b,r,param)
          *BasisPoly(p-1,r,-2.0*xx/b+(b-1.0)/b);
    }
    return sum;
  }
  else if ((xx>-a) && (xx<a)) {
    for (r=0; r<p; r++) {
      sum+=pow(a,(double)r)
          *( k(-a,r,param)*BasisPoly(p-1,r,xx/a)
              +k( a,r,param)*BasisPoly(p-1,r,-xx/a)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx>0.5-b) {
    sum=(k(-0.5,0,param)+k(0.5,0,param))/2.0
        *BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b);
    for (r=0; r<p; r++) {
      sum+=pow(b/2.0,(double)r)
          *k(0.5-b,r,param)
          *BasisPoly(p-1,r,2.0*xx/b-(1.0-b)/b);
    }
    return sum;
  }
  return k(xx,0,param);
}

/** regularized kernel with K_I arbitrary and K_B periodized
 *  (used in 1D)
 */
double _Complex regkern1(kernel k, double xx, int p, const double *param, double a, double b)
{
  int r;
  double _Complex sum=0.0;

  if (xx<-0.5)
    xx=-0.5;
  if (xx>0.5)
    xx=0.5;
  if ((xx>=-0.5+b && xx<=-a) || (xx>=a && xx<=0.5-b))
  {
    return k(xx,0,param);
  }
  else if ((xx>-a) && (xx<a))
  {
    for (r=0; r<p; r++) {
      sum+=pow(a,(double)r)
          *( k(-a,r,param)*BasisPoly(p-1,r,xx/a)
              +k( a,r,param)*BasisPoly(p-1,r,-xx/a)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx<-0.5+b)
  {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)
          *( k(0.5-b,r,param)*BasisPoly(p-1,r,(xx+0.5)/b)
              +k(-0.5+b,r,param)*BasisPoly(p-1,r,-(xx+0.5)/b)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  else if (xx>0.5-b)
  {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)
          *( k(0.5-b,r,param)*BasisPoly(p-1,r,(xx-0.5)/b)
              +k(-0.5+b,r,param)*BasisPoly(p-1,r,-(xx-0.5)/b)*(r & 1 ? -1 : 1));
    }
    return sum;
  }
  return k(xx,0,param);
}

/** regularized kernel for even kernels with K_I even and K_B mirrored */
double _Complex regkern2(kernel k, double xx, int p, const double *param, double a, double b)
{
  int r;
  double _Complex sum=0.0;

  xx=fabs(xx);

  if (xx>0.5) {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)*k(0.5-b,r,param)
          *(BasisPoly(p-1,r,0)+BasisPoly(p-1,r,0));
    }
    return sum;
  }
  else if ((a<=xx) && (xx<=0.5-b)) {
    return k(xx,0,param);
  }
  else if (xx<a) {
    for (r=0; r<p; r++) {
      sum+=pow(-a,(double)r)*k(a,r,param)
          *(BasisPoly(p-1,r,xx/a)+BasisPoly(p-1,r,-xx/a));
    }
    return sum;
  }
  else if ((0.5-b<xx) && (xx<=0.5)) {
    for (r=0; r<p; r++) {
      sum+=pow(b,(double)r)*k(0.5-b,r,param)
          *(BasisPoly(p-1,r,(xx-0.5)/b)+BasisPoly(p-1,r,-(xx-0.5)/b));
    }
    return sum;
  }
  return 0.0;
}

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth to K(1/2) (used in dD, d>1)
 */
double _Complex regkern3(kernel k, double xx, int p, const double *param, double a, double b)
{
  int r;
  double _Complex sum=0.0;

  xx=fabs(xx);

  if (xx>=0.5) {
    /*return kern(typ,c,0,0.5);*/
    xx=0.5;
  }
  /* else */
  if ((a<=xx) && (xx<=0.5-b)) {
    return k(xx,0,param);
  }
  else if (xx<a) {
    for (r=0; r<p; r++) {
      sum+=pow(-a,(double)r)*k(a,r,param)
          *(BasisPoly(p-1,r,xx/a)+BasisPoly(p-1,r,-xx/a));
    }
    /*sum=kern(typ,c,0,xx); */
    return sum;
  }
  else if ((0.5-b<xx) && (xx<=0.5)) {
    sum=k(0.5,0,param)*BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b);
    /* sum=regkern2(typ,c,p,a,b, 0.5)*BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b); */
    for (r=0; r<p; r++) {
      sum+=pow(b/2.0,(double)r)
          *k(0.5-b,r,param)
          *BasisPoly(p-1,r,2.0*xx/b-(1.0-b)/b);
    }
    return sum;
  }
  return 0.0;
}

/** linear spline interpolation in near field with even kernels */
double _Complex linintkern(const double x, const double _Complex *Add,
  const int Ad, const double a)
{
  double c,c1,c3;
  int r;
  double _Complex f1,f2;

  c=x*Ad/a;
  r=c; r=abs(r);
  f1=Add[r];f2=Add[r+1];
  c=fabs(c);
  c1=c-r;
  c3=c1-1.0;
  return (-f1*c3+f2*c1);
}

double _Complex quadrintkern(const double x, const double _Complex *Add,
  const int Ad, const double a)
{
  double c,c1,c2,c3;
  int r;
  double _Complex f0,f1,f2;

  c=x*Ad/a;
  r=c; r=abs(r);
  if (r==0) {f0=Add[r+1];f1=Add[r];f2=Add[r+1];}
  else { f0=Add[r-1];f1=Add[r];f2=Add[r+1];}
  c=fabs(c);
  c1=c-r;
  c2=c1+1.0;
  c3=c1-1.0;
  return (f0*c1*c3/2.0-f1*c2*c3+f2*c2*c1/2.0);
}

/** cubic spline interpolation in near field with even kernels */
double _Complex kubintkern(const double x, const double _Complex *Add,
  const int Ad, const double a)
{
  double c,c1,c2,c3,c4;
  int r;
  double _Complex f0,f1,f2,f3;
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
double _Complex kubintkern1(const double x, const double _Complex *Add,
  const int Ad, const double a)
{
  double c,c1,c2,c3,c4;
  int r;
  double _Complex f0,f1,f2,f3;
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

#ifdef NF_BO
/** initialize box-based search data structures */
void BuildBox(fastsum_plan *ths)
{
  int t, l;
  int *box_index;
  double val[ths->d];

  box_index = (int *) malloc(ths->box_count * sizeof(int));
  for (t=0; t < ths->box_count; t++)
    box_index[t] = 0;

  for (l=0; l < ths->N_total; l++)
  {
    int ind = 0;
    for (t=0; t < ths->d; t++)
    {
      val[t] = ths->x[ths->d * l + t] + 0.25 - ths->eps_B/2.0;
      ind *= ths->box_count_per_dim;
      ind += (int) (val[t] / ths->eps_I);
    }
    box_index[ind]++;
  }

  ths->box_offset[0] = 0;
  for (t=1; t<=ths->box_count; t++)
  {
    ths->box_offset[t] = ths->box_offset[t-1] + box_index[t-1];
    box_index[t-1] = ths->box_offset[t-1];
  }

  for (l=0; l < ths->N_total; l++)
  {
    int ind = 0;
    for (t=0; t < ths->d; t++)
    {
      val[t] = ths->x[ths->d * l + t] + 0.25 - ths->eps_B/2.0;
      ind *= ths->box_count_per_dim;
      ind += (int) (val[t] / ths->eps_I);
    }

    ths->box_alpha[box_index[ind]] = ths->alpha[l];

    for (t=0; t < ths->d; t++)
    {
      ths->box_x[ths->d * box_index[ind] + t] = ths->x[ths->d * l + t];
    }
    box_index[ind]++;
  }
  free(box_index);
}

/** inner computation function for box-based near field correction */
inline double _Complex calc_SearchBox(int d, double *y, double *x, double _Complex *alpha, int start, int end_lt, const double _Complex *Add, const int Ad, int p, double a, const kernel k, const double *param, const unsigned flags)
{
  double _Complex result = 0.0;

  int m, l;
  double r;

  for (m = start; m < end_lt; m++)
  {
      if (d==1)
      {
        r = y[0]-x[m];
      }
      else
      {
        r=0.0;
        for (l=0; l<d; l++)
          r+=(y[l]-x[m*d+l])*(y[l]-x[m*d+l]);
        r=sqrt(r);
      }
      if (fabs(r)<a)
      {
        result += alpha[m]*k(r,0,param);              /* alpha*(kern-regkern) */
	if (d==1)
	{
          if (flags & EXACT_NEARFIELD)
            result -= alpha[m]*regkern1(k,r,p,param,a,1.0/16.0); /* exact value (in 1D)  */
          else
            result -= alpha[m]*kubintkern1(r,Add,Ad,a);               /* spline approximation */
	}
	else
	{
          if (flags & EXACT_NEARFIELD)
            result -= alpha[m]*regkern(k,r,p,param,a,1.0/16.0);  /* exact value (in dD)  */
          else
#if defined(NF_KUB)
            result -= alpha[m]*kubintkern(r,Add,Ad,a);                /* spline approximation */
#elif defined(NF_QUADR)
            result -= alpha[m]*quadrintkern(r,Add,Ad,a);                /* spline approximation */
#elif defined(NF_LIN)
            result -= alpha[m]*linintkern(r,Add,Ad,a);                /* spline approximation */
#else
  #error define interpolation method
#endif
        }
      }
  }
  return result;
}

/** box-based near field correction */
double _Complex SearchBox(double *y, fastsum_plan *ths)
{
  double _Complex val = 0.0;
  int t, l;
  int y_multiind[ths->d];
  int multiindex[ths->d];
  int y_ind;

  for (t=0; t < ths->d; t++)
  {
    y_multiind[t] = ((y[t] + 0.25 - ths->eps_B/2.0) / ths->eps_I);
  }

  if (ths->d==1)
  {
    for (y_ind = max_i(0, y_multiind[0]-1); y_ind < ths->box_count_per_dim && y_ind <= y_multiind[0]+1; y_ind++){
      val += calc_SearchBox(ths->d, y, ths->box_x, ths->box_alpha, ths->box_offset[y_ind], ths->box_offset[y_ind+1], ths->Add, ths->Ad, ths->p, ths->eps_I, ths->k, ths->kernel_param, ths->flags);
      }
  }
  else if (ths->d==2)
  {
    for (multiindex[0] = max_i(0, y_multiind[0]-1); multiindex[0] < ths->box_count_per_dim && multiindex[0] <= y_multiind[0]+1; multiindex[0]++)
      for (multiindex[1] = max_i(0, y_multiind[1]-1); multiindex[1] < ths->box_count_per_dim && multiindex[1] <= y_multiind[1]+1; multiindex[1]++)
      {
        y_ind = (ths->box_count_per_dim * multiindex[0]) + multiindex[1];
        val += calc_SearchBox(ths->d, y, ths->box_x, ths->box_alpha, ths->box_offset[y_ind], ths->box_offset[y_ind+1], ths->Add, ths->Ad, ths->p, ths->eps_I, ths->k, ths->kernel_param, ths->flags);
      }
  }
  else if(ths->d==3)
  {
    for (multiindex[0] = max_i(0, y_multiind[0]-1); multiindex[0] < ths->box_count_per_dim && multiindex[0] <= y_multiind[0]+1; multiindex[0]++)
      for (multiindex[1] = max_i(0, y_multiind[1]-1); multiindex[1] < ths->box_count_per_dim && multiindex[1] <= y_multiind[1]+1; multiindex[1]++)
        for (multiindex[2] = max_i(0, y_multiind[2]-1); multiindex[2] < ths->box_count_per_dim && multiindex[2] <= y_multiind[2]+1; multiindex[2]++)
        {
          y_ind = ((ths->box_count_per_dim * multiindex[0]) + multiindex[1]) * ths->box_count_per_dim + multiindex[2];
          val += calc_SearchBox(ths->d, y, ths->box_x, ths->box_alpha, ths->box_offset[y_ind], ths->box_offset[y_ind+1], ths->Add, ths->Ad, ths->p, ths->eps_I, ths->k, ths->kernel_param, ths->flags);
        }
  }
  else {
    exit(-1);
  }
  return val;
}
#endif

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
double _Complex SearchTree(const int d, const int t, const double *x,
  const double _Complex *alpha, const double *xmin, const double *xmax,
  const int N, const kernel k, const double *param, const int Ad,
  const double _Complex *Add, const int p, const unsigned flags)
{
  int m=N/2;
  double Min=xmin[t], Max=xmax[t], Median=x[m*d+t];
  double a=fabs(Max-Min)/2;
  int l;
  int E=0;
  double r;

  if (N==0)
    return 0.0;
  if (Min>Median)
    return SearchTree(d,(t+1)%d,x+(m+1)*d,alpha+(m+1),xmin,xmax,N-m-1,k,param,Ad,Add,p,flags);
  else if (Max<Median)
    return SearchTree(d,(t+1)%d,x,alpha,xmin,xmax,m,k,param,Ad,Add,p,flags);
  else
  {
    double _Complex result = 0.0;
    E=0;

    for (l=0; l<d; l++)
    {
      if ( x[m*d+l]>xmin[l] && x[m*d+l]<xmax[l] )
        E++;
    }

    if (E==d)
    {
      if (d==1)
      {
        r = xmin[0]+a-x[m];  /* remember: xmin+a = y */
      }
      else
      {
        r=0.0;
        for (l=0; l<d; l++)
          r+=(xmin[l]+a-x[m*d+l])*(xmin[l]+a-x[m*d+l]);  /* remember: xmin+a = y */
        r=sqrt(r);
      }
      if (fabs(r)<a)
      {
        result += alpha[m]*k(r,0,param);                         /* alpha*(kern-regkern) */
        if (d==1)
        {
          if (flags & EXACT_NEARFIELD)
            result -= alpha[m]*regkern1(k,r,p,param,a,1.0/16.0); /* exact value (in 1D)  */
          else
            result -= alpha[m]*kubintkern1(r,Add,Ad,a);               /* spline approximation */
        }
        else
        {
          if (flags & EXACT_NEARFIELD)
            result -= alpha[m]*regkern(k,r,p,param,a,1.0/16.0);  /* exact value (in dD)  */
          else
#if defined(NF_KUB)
            result -= alpha[m]*kubintkern(r,Add,Ad,a);                /* spline approximation */
#elif defined(NF_QUADR)
            result -= alpha[m]*quadrintkern(r,Add,Ad,a);                /* spline approximation */
#elif defined(NF_LIN)
            result -= alpha[m]*linintkern(r,Add,Ad,a);                /* spline approximation */
#else
  #error define interpolation method
#endif
        }
      }
    }
    result += SearchTree(d,(t+1)%d,x+(m+1)*d,alpha+(m+1),xmin,xmax,N-m-1,k,param,Ad,Add,p,flags)
      + SearchTree(d,(t+1)%d,x,alpha,xmin,xmax,m,k,param,Ad,Add,p,flags);
    return result;
  }
}

/** initialization of fastsum plan */
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total, kernel k, double *param, unsigned flags, int nn, int m, int p, double eps_I, double eps_B)
{
  int t;
  int N[d], n[d];
  int n_total;
  int sort_flags_trafo = 0;
  int sort_flags_adjoint = 0;
#ifdef _OPENMP
  int nthreads = nfft_get_omp_num_threads();
#endif

  if (d > 1)
  {
    sort_flags_trafo = NFFT_SORT_NODES;
#ifdef _OPENMP
    sort_flags_adjoint = NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT;
#else
    sort_flags_adjoint = NFFT_SORT_NODES;
#endif
  }

  ths->d = d;

  ths->N_total = N_total;
  ths->M_total = M_total;

  ths->x = (double *)nfft_malloc(d*N_total*(sizeof(double)));
  ths->alpha = (double _Complex *)nfft_malloc(N_total*(sizeof(double _Complex)));

  ths->y = (double *)nfft_malloc(d*M_total*(sizeof(double)));
  ths->f = (double _Complex *)nfft_malloc(M_total*(sizeof(double _Complex)));

  ths->k = k;
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
      ths->Add = (double _Complex *)nfft_malloc((ths->Ad+5)*(sizeof(double _Complex)));
    }
    else
    {
      if (ths->k == one_over_x)
      {
        double delta = 1e-8;
        switch(p)
        {
          case 2: delta = 1e-3;
                  break;
          case 3: delta = 1e-4;
                  break;
          case 4: delta = 1e-5;
                  break;
          case 5: delta = 1e-6;
                  break;
          case 6: delta = 1e-6;
                  break;
          case 7: delta = 1e-7;
                  break;
          default: delta = 1e-8;
        }

#if defined(NF_KUB)
        ths->Ad = max_i(10, (int) ceil(1.4/pow(delta,1.0/4.0)));
        ths->Add = (double _Complex *)nfft_malloc((ths->Ad+3)*(sizeof(double _Complex)));
#elif defined(NF_QUADR)
        ths->Ad = (int) ceil(2.2/pow(delta,1.0/3.0));
        ths->Add = (double _Complex *)nfft_malloc((ths->Ad+3)*(sizeof(double _Complex)));
#elif defined(NF_LIN)
        ths->Ad = (int) ceil(1.7/pow(delta,1.0/2.0));
        ths->Add = (double _Complex *)nfft_malloc((ths->Ad+3)*(sizeof(double _Complex)));
#else
  #error define NF_LIN or NF_QUADR or NF_KUB
#endif
      }
      else
      {
        ths->Ad = 2*(ths->p)*(ths->p);
        ths->Add = (double _Complex *)nfft_malloc((ths->Ad+3)*(sizeof(double _Complex)));
      }
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
                   sort_flags_adjoint |
                   PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                   FFTW_MEASURE| FFTW_DESTROY_INPUT);
  nfft_init_guru(&(ths->mv2), d, N, M_total, n, m,
                   sort_flags_trafo |
                   PRE_PHI_HUT| PRE_PSI| MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                   FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /** init d-dimensional FFTW plan */
  n_total = 1;
  for (t=0; t<d; t++)
    n_total *= nn;

  ths->b = (fftw_complex *)nfft_malloc(n_total*sizeof(fftw_complex));
#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
  fftw_plan_with_nthreads(nthreads);
#endif

  ths->fft_plan = fftw_plan_dft(d,N,ths->b,ths->b,FFTW_FORWARD,FFTW_ESTIMATE);

#ifdef _OPENMP
}
#endif

#ifdef NF_BO
  ths->box_count_per_dim = floor((0.5 - ths->eps_B) / ths->eps_I) + 1;
  ths->box_count = 1;
  for (t=0; t<ths->d; t++)
    ths->box_count *= ths->box_count_per_dim;

  ths->box_offset = (int *) nfft_malloc((ths->box_count+1) * sizeof(int));

  ths->box_alpha = (double _Complex *)nfft_malloc(ths->N_total*(sizeof(double _Complex)));

  ths->box_x = (double *) nfft_malloc(ths->d * ths->N_total *  sizeof(double));
#endif
}

/** finalization of fastsum plan */
void fastsum_finalize(fastsum_plan *ths)
{
  nfft_free(ths->x);
  nfft_free(ths->alpha);
  nfft_free(ths->y);
  nfft_free(ths->f);

  if (!(ths->flags & EXACT_NEARFIELD))
    nfft_free(ths->Add);

  nfft_finalize(&(ths->mv1));
  nfft_finalize(&(ths->mv2));

#ifdef _OPENMP
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
#endif
  fftw_destroy_plan(ths->fft_plan);
#ifdef _OPENMP
}
#endif

  nfft_free(ths->b);

#ifdef NF_BO
  nfft_free(ths->box_offset);
  nfft_free(ths->box_alpha);
  nfft_free(ths->box_x);
#endif
}

/** direct computation of sums */
void fastsum_exact(fastsum_plan *ths)
{
  int j,k;
  int t;
  double r;

  #pragma omp parallel for default(shared) private(j,k,t,r)
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
      ths->f[j] += ths->alpha[k] * ths->k(r,0,ths->kernel_param);
    }
  }
}

/** precomputation for fastsum */
void fastsum_precompute(fastsum_plan *ths)
{
  int j,k,t;
  int n_total;
  ticks t0, t1;

  ths->MEASURE_TIME_t[0] = 0.0;
  ths->MEASURE_TIME_t[1] = 0.0;
  ths->MEASURE_TIME_t[2] = 0.0;
  ths->MEASURE_TIME_t[3] = 0.0;

#ifdef MEASURE_TIME
  t0 = getticks();
#endif

#if defined(NF_ST)
  /** sort source knots */
  BuildTree(ths->d,0,ths->x,ths->alpha,ths->N_total);
#elif defined(NF_BO)  
  BuildBox(ths);
#else
  #error Either define NF_ST or NF_BO
#endif

#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[3] += nfft_elapsed_seconds(t1,t0);
#endif


#ifdef MEASURE_TIME
  t0 = getticks();
#endif
  /** precompute spline values for near field*/
  if (!(ths->flags & EXACT_NEARFIELD))
  {
    if (ths->d==1)
      #pragma omp parallel for default(shared) private(k)
      for (k=-ths->Ad/2-2; k <= ths->Ad/2+2; k++)
        ths->Add[k+ths->Ad/2+2] = regkern1(ths->k, ths->eps_I*(double)k/ths->Ad*2, ths->p, ths->kernel_param, ths->eps_I, ths->eps_B);
    else
      #pragma omp parallel for default(shared) private(k)
      for (k=0; k <= ths->Ad+2; k++)
        ths->Add[k] = regkern3(ths->k, ths->eps_I*(double)k/ths->Ad, ths->p, ths->kernel_param, ths->eps_I, ths->eps_B);
  }
#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[0] += nfft_elapsed_seconds(t1,t0);
#endif


#ifdef MEASURE_TIME
  t0 = getticks();
#endif
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
#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[1] += nfft_elapsed_seconds(t1,t0);
#endif

  /** init Fourier coefficients */
  for(k=0; k<ths->mv1.M_total;k++)
    ths->mv1.f[k] = ths->alpha[k];

#ifdef MEASURE_TIME
  t0 = getticks();
#endif
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
#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[2] += nfft_elapsed_seconds(t1,t0);
#endif


#ifdef MEASURE_TIME
  t0 = getticks();
#endif
  /** precompute Fourier coefficients of regularised kernel*/
  n_total = 1;
  for (t=0; t<ths->d; t++)
    n_total *= ths->n;

  #pragma omp parallel for default(shared) private(j,k,t)
  for (j=0; j<n_total; j++)
  {
    if (ths->d==1)
      ths->b[j] = regkern1(ths->k, (double)j / (ths->n) - 0.5, ths->p, ths->kernel_param, ths->eps_I, ths->eps_B)/n_total;
    else
    {
      k=j;
      ths->b[j]=0.0;
      for (t=0; t<ths->d; t++)
      {
        ths->b[j] += ((double)(k % (ths->n)) / (ths->n) - 0.5) * ((double)(k % (ths->n)) / (ths->n) - 0.5);
        k = k / (ths->n);
      }
      ths->b[j] = regkern3(ths->k, sqrt(ths->b[j]), ths->p, ths->kernel_param, ths->eps_I, ths->eps_B)/n_total;
    }
  }

  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);
  fftw_execute(ths->fft_plan);
  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);
#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[0] += nfft_elapsed_seconds(t1,t0);
#endif
}

/** fast NFFT-based summation */
void fastsum_trafo(fastsum_plan *ths)
{
  int j,k,t;
  ticks t0, t1;

  ths->MEASURE_TIME_t[4] = 0.0; 
  ths->MEASURE_TIME_t[5] = 0.0;
  ths->MEASURE_TIME_t[6] = 0.0;
  ths->MEASURE_TIME_t[7] = 0.0;

#ifdef MEASURE_TIME
  t0 = getticks();
#endif
  /** first step of algorithm */
  nfft_adjoint(&(ths->mv1));
#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[4] += nfft_elapsed_seconds(t1,t0);
#endif


#ifdef MEASURE_TIME
  t0 = getticks();
#endif
  /** second step of algorithm */
  #pragma omp parallel for default(shared) private(k)
  for (k=0; k<ths->mv2.N_total; k++)
    ths->mv2.f_hat[k] = ths->b[k] * ths->mv1.f_hat[k];
#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[5] += nfft_elapsed_seconds(t1,t0);
#endif


#ifdef MEASURE_TIME
  t0 = getticks();
#endif
  /** third step of algorithm */
  nfft_trafo(&(ths->mv2));
#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[6] += nfft_elapsed_seconds(t1,t0);
#endif


#ifdef MEASURE_TIME
  t0 = getticks();
#endif
  /** add near field */
  #pragma omp parallel for default(shared) private(j,k,t)
  for (j=0; j<ths->M_total; j++)
  {
    double ymin[ths->d], ymax[ths->d]; /** limits for d-dimensional near field box */
#if defined(NF_ST)
    for (t=0; t<ths->d; t++)
    {
      ymin[t] = ths->y[ths->d*j+t] - ths->eps_I;
      ymax[t] = ths->y[ths->d*j+t] + ths->eps_I;
    }
    ths->f[j] = ths->mv2.f[j] + SearchTree(ths->d,0, ths->x, ths->alpha, ymin, ymax, ths->N_total, ths->k, ths->kernel_param, ths->Ad, ths->Add, ths->p, ths->flags);
#elif defined(NF_BO)
    ths->f[j] = ths->mv2.f[j] + SearchBox(ths->y + ths->d*j, ths);
#else
  #error missing NF_ST or NF_BO
#endif
    /* ths->f[j] = ths->mv2.f[j]; */
    /* ths->f[j] = SearchTree(ths->d,0, ths->x, ths->alpha, ymin, ymax, ths->N_total, ths->k, ths->kernel_param, ths->Ad, ths->Add, ths->p, ths->flags); */
  }

#ifdef MEASURE_TIME
  t1 = getticks();
  ths->MEASURE_TIME_t[7] += nfft_elapsed_seconds(t1,t0);
#endif
}
/* \} */

/* fastsum.c */
