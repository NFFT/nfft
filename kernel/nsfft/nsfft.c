/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include "nfft3.h"
#include "infft.h"

#define NSFTT_DISABLE_TEST

/* computes a 2d ndft by 1d nfft along the dimension 1 times
   1d ndft along dimension 0
*/
static void short_nfft_trafo_2d(nfft_plan* ths, nfft_plan* plan_1d)
{
  int j,k0;
  double omega;

  for(j=0;j<ths->M_total;j++)
    {
      ths->f[j]= 0;
      plan_1d->x[j] = ths->x[ths->d * j + 1];
    }

  for(k0=0;k0<ths->N[0];k0++) /* for shorties */
    {
      plan_1d->f_hat = ths->f_hat + k0*ths->N[1];

      nfft_trafo(plan_1d);

      for(j=0;j<ths->M_total;j++)
	{
	  omega = ((double)(k0 - ths->N[0]/2)) * ths->x[ths->d * j + 0];
          ths->f[j] += plan_1d->f[j] * cexp( - I*2*KPI*omega);
	}
    }
}

static void short_nfft_adjoint_2d(nfft_plan* ths, nfft_plan* plan_1d)
{
  int j,k0;
  double omega;

  for(j=0;j<ths->M_total;j++)
    plan_1d->x[j] = ths->x[ths->d * j + 1];

  for(k0=0;k0<ths->N[0];k0++) /* for shorties */
    {
      for(j=0;j<ths->M_total;j++)
	{
	  omega = ((double)(k0 - ths->N[0]/2)) * ths->x[ths->d * j + 0];
          plan_1d->f[j] = ths->f[j] * cexp( + _Complex_I*2*KPI*omega);
	}

      plan_1d->f_hat = ths->f_hat + k0*ths->N[1];

      nfft_adjoint(plan_1d);
    }
}

/* computes a 3d ndft by 1d nfft along the dimension 2 times
   2d ndft along dimension 0,1
*/
static void short_nfft_trafo_3d_1(nfft_plan* ths, nfft_plan* plan_1d)
{
  int j,k0,k1;
  double omega;

  for(j=0;j<ths->M_total;j++)
    {
      ths->f[j] = 0;
      plan_1d->x[j] = ths->x[ths->d * j + 2];
    }

  for(k0=0;k0<ths->N[0];k0++) /* for shorties */
    for(k1=0;k1<ths->N[1];k1++)
      {
	plan_1d->f_hat = ths->f_hat + (k0*ths->N[1]+k1)*ths->N[2];

	nfft_trafo(plan_1d);

	for(j=0;j<ths->M_total;j++)
	  {
	    omega = ((double)(k0 - ths->N[0]/2)) * ths->x[ths->d * j + 0]
	      +     ((double)(k1 - ths->N[1]/2)) * ths->x[ths->d * j + 1];
            ths->f[j] += plan_1d->f[j] * cexp( - I*2*KPI*omega);
	  }
      }
}

static void short_nfft_adjoint_3d_1(nfft_plan* ths, nfft_plan* plan_1d)
{
  int j,k0,k1;
  double omega;

  for(j=0;j<ths->M_total;j++)
    plan_1d->x[j] = ths->x[ths->d * j + 2];

  for(k0=0;k0<ths->N[0];k0++) /* for shorties */
    for(k1=0;k1<ths->N[1];k1++)
      {
	for(j=0;j<ths->M_total;j++)
	  {
	    omega = ((double)(k0 - ths->N[0]/2)) * ths->x[ths->d * j + 0]
	      +     ((double)(k1 - ths->N[1]/2)) * ths->x[ths->d * j + 1];
            plan_1d->f[j] = ths->f[j] * cexp( + _Complex_I*2*KPI*omega);
	  }

	plan_1d->f_hat = ths->f_hat + (k0*ths->N[1]+k1)*ths->N[2];

	nfft_adjoint(plan_1d);
      }
}

/* computes a 3d ndft by 2d nfft along the dimension 1,2 times
   1d ndft along dimension 0
*/
static void short_nfft_trafo_3d_2(nfft_plan* ths, nfft_plan* plan_2d)
{
  int j,k0;
  double omega;

  for(j=0;j<ths->M_total;j++)
    {
      ths->f[j] = 0;
      plan_2d->x[2*j+0] = ths->x[ths->d * j + 1];
      plan_2d->x[2*j+1] = ths->x[ths->d * j + 2];
    }

  for(k0=0;k0<ths->N[0];k0++) /* for shorties */
    {
      plan_2d->f_hat = ths->f_hat + k0*ths->N[1]*ths->N[2];

      nfft_trafo(plan_2d);

      for(j=0;j<ths->M_total;j++)
	{
	  omega = ((double)(k0 - ths->N[0]/2)) * ths->x[ths->d * j + 0];
	  ths->f[j] += plan_2d->f[j] * cexp( - I*2*KPI*omega);
	}
    }
}

static void short_nfft_adjoint_3d_2(nfft_plan* ths, nfft_plan* plan_2d)
{
  int j,k0;
  double omega;

  for(j=0;j<ths->M_total;j++)
    {
      plan_2d->x[2*j+0] = ths->x[ths->d * j + 1];
      plan_2d->x[2*j+1] = ths->x[ths->d * j + 2];
    }

  for(k0=0;k0<ths->N[0];k0++) /* for shorties */
    {
      for(j=0;j<ths->M_total;j++)
	{
	  omega = ((double)(k0 - ths->N[0]/2)) * ths->x[ths->d * j + 0];
	  plan_2d->f[j] = ths->f[j] * cexp( + _Complex_I*2*KPI*omega);
	}

      plan_2d->f_hat = ths->f_hat + k0*ths->N[1]*ths->N[2];

      nfft_adjoint(plan_2d);
    }
}

/*---------------------------------------------------------------------------*/

#ifdef GAUSSIAN
static int index_sparse_to_full_direct_2d(int J, int k)
{
    int N=X(exp2i)(J+2);               /* number of full coeffs             */
    int N_B=X(exp2i)(J);               /* number in each sparse block       */

    int j=k/N_B;                        /* consecutive number of Block       */
    int r=j/4;                          /* level of block                    */

    int i, o, a, b,s,l,m1,m2;
    int k1,k2;

    if (k>=(J+4)*X(exp2i)(J+1))
      {
	printf("Fehler!\n");
	return(-1);
      }
    else
      {
	if (r>(J+1)/2)                  /* center block                      */
	  {
	    i=k-4*((J+1)/2+1)*N_B;
	    a=X(exp2i)(J/2+1);
	    m1=i/a;
	    m2=i%a;
	    k1=N/2-a/2+m1;
	    k2=N/2-a/2+m2;
	  }
	else                            /* no center block                   */
	  {
	    i=k-j*N_B;                  /* index in specific block           */
	    o=j%4;                      /* kind of specific block            */
	    a=X(exp2i)(r);
	    b=X(exp2i)(J-r);
	    l=MAX(a,b);                 /* long dimension of block           */
	    s=MIN(a,b);                 /* short dimension of block          */
	    m1=i/l;
	    m2=i%l;

	    switch(o)
	      {
	      case 0:
		{
		  k1=N/2-a/2 ;
		  k2=N/2+ b  ;

		  if (b>=a)
		    {
		      k1+=m1;
		      k2+=m2;
		    }
		  else
		    {
		      k1+=m2;
		      k2+=m1;
		    }
		  break;
		}
	      case 1:
		{
		  k1=N/2+ b  ;
		  k2=N/2-a/2 ;

		  if (b>a)
		    {
		      k1+=m2;
		      k2+=m1;
		    }
		  else
		    {
		      k1+=m1;
		      k2+=m2;
		    }
		  break;
		}
	      case 2:
		{
		  k1=N/2-a/2 ;
		  k2=N/2-2*b ;

		  if (b>=a)
		    {
		      k1+=m1;
		      k2+=m2;
		    }
		  else
		    {
		      k1+=m2;
		      k2+=m1;
		    }
		  break;
		}
	      case 3:
		{
		  k1=N/2-2*b ;
		  k2=N/2-a/2 ;

		  if (b>a)
		    {
		      k1+=m2;
		      k2+=m1;
		    }
		  else
		    {
		      k1+=m1;
		      k2+=m2;
		    }
		  break;
		}
	      default:
		{
		  k1=-1;
		  k2=-1;
		}
	      }
	  }
	//printf("m1=%d, m2=%d\n",m1,m2);
	return(k1*N+k2);
      }
}
#endif

static inline int index_sparse_to_full_2d(nsfft_plan *ths, int k)
{
  /* only by lookup table */
  if( k < ths->N_total)
    return ths->index_sparse_to_full[k];
  else
    return -1;
}

#ifndef NSFTT_DISABLE_TEST
static int index_full_to_sparse_2d(int J, int k)
{
    int N=X(exp2i)(J+2);               /* number of full coeffs       */
    int N_B=X(exp2i)(J);               /* number in each sparse block */

    int k1=k/N-N/2;                     /* coordinates in the full grid */
    int k2=k%N-N/2;                     /* k1: row, k2: column          */

    int r,a,b;

    a=X(exp2i)(J/2+1);

    if ( (k1>=-(a/2)) && (k1<a/2) && (k2>=(-a/2)) && (k2<a/2) )
      {
	return(4*((J+1)/2+1)*N_B+(k1+a/2)*a+(k2+a/2));
      }

    for (r=0; r<=(J+1)/2; r++)
      {
	b=X(exp2i)(r);
	a=X(exp2i)(J-r);
	if ( (k1>=-(b/2)) && (k1<(b+1)/2) && (k2>=a) && (k2<2*a) )
	  {
            if (a>=b)
	      return((4*r+0)*N_B+(k1+b/2)*a+(k2-a));
	    else
	      return((4*r+0)*N_B+(k2-a)*b+(k1+b/2));
	  }
	else if ( (k1>=a) && (k1<2*a) && (k2>=-(b/2)) && (k2<(b+1)/2) )
	  {
            if (a>b)
	      return((4*r+1)*N_B+(k2+b/2)*a+(k1-a));
	    else
	      return((4*r+1)*N_B+(k1-a)*b+(k2+b/2));
	  }
	else if ( (k1>=-(b/2)) && (k1<(b+1)/2) && (k2>=-2*a) && (k2<-a) )
	  {
            if (a>=b)
	      return((4*r+2)*N_B+(k1+b/2)*a+(k2+2*a));
	    else
	      return((4*r+2)*N_B+(k2+2*a)*b+(k1+b/2));
	  }
	else if ( (k1>=-2*a) && (k1<-a) && (k2>=-(b/2)) && (k2<(b+1)/2) )
	  {
            if (a>b)
	      return((4*r+3)*N_B+(k2+b/2)*a+(k1+2*a));
	    else
	      return((4*r+3)*N_B+(k1+2*a)*b+(k2+b/2));
	  }
      }

    return(-1);
}
#endif

#ifdef GAUSSIAN
static void init_index_sparse_to_full_2d(nsfft_plan *ths)
{
  int k_S;

  for (k_S=0; k_S<ths->N_total; k_S++)
    ths->index_sparse_to_full[k_S]=index_sparse_to_full_direct_2d(ths->J, k_S);
}
#endif

#ifdef GAUSSIAN
static inline int index_sparse_to_full_3d(nsfft_plan *ths, int k)
{
  /* only by lookup table */
  if( k < ths->N_total)
    return ths->index_sparse_to_full[k];
  else
    return -1;
}
#endif

#ifndef NSFTT_DISABLE_TEST
static int index_full_to_sparse_3d(int J, int k)
{
  int N=X(exp2i)(J+2);                 /* length of the full grid           */
  int N_B_r;                            /* size of a sparse block in level r */
  int sum_N_B_less_r;                   /* sum N_B_r                         */

  int r,a,b;

  int k3=(k%N)-N/2;                     /* coordinates in the full grid      */
  int k2=((k/N)%N)-N/2;
  int k1=k/(N*N)-N/2;

  a=X(exp2i)(J/2+1);                   /* length of center block            */

  if((k1>=-(a/2)) && (k1<a/2) && (k2>=(-a/2)) && (k2<a/2) && (k3>=(-a/2)) &&
     (k3<a/2))
    {
      return(6*X(exp2i)(J)*(X(exp2i)((J+1)/2+1)-1)+((k1+a/2)*a+(k2+a/2))*a+
             (k3+a/2));
    }

  sum_N_B_less_r=0;
  for (r=0; r<=(J+1)/2; r++)
    {
      a=X(exp2i)(J-r);
      b=X(exp2i)(r);

      N_B_r=a*b*b;

      /* right - rear - top - left - front - bottom */
      if ((k1>=a) && (k1<2*a) && (k2>=-(b/2)) && (k2<(b+1)/2) &&
          (k3>=-(b/2)) && (k3<(b+1)/2)) /* right */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*0 + ((k2+b/2)*b+k3+b/2)*a + (k1-a);
	  else
	    return sum_N_B_less_r+N_B_r*0 + ((k1-a)*b+(k2+b/2))*b + (k3+b/2);
	}
      else if ((k2>=a) && (k2<2*a) && (k1>=-(b/2)) && (k1<(b+1)/2) &&
               (k3>=-(b/2)) && (k3<(b+1)/2)) /* rear */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*1 + ((k1+b/2)*b+k3+b/2)*a + (k2-a);
	  else if (a==b)
	    return sum_N_B_less_r+N_B_r*1 + ((k1+b/2)*b+(k2-a))*a + (k3+b/2);
	  else
	    return sum_N_B_less_r+N_B_r*1 + ((k2-a)*b+(k1+b/2))*b + (k3+b/2);
	}
       else if ((k3>=a) && (k3<2*a) && (k1>=-(b/2)) && (k1<(b+1)/2) &&
                (k2>=-(b/2)) && (k2<(b+1)/2)) /* top */
	{
	  if(a>=b)
	    return sum_N_B_less_r+N_B_r*2 + ((k1+b/2)*b+k2+b/2)*a + (k3-a);
	  else
	    return sum_N_B_less_r+N_B_r*2 + ((k3-a)*b+(k1+b/2))*b + (k2+b/2);
	}

      else if ((k1>=-2*a) && (k1<-a) && (k2>=-(b/2)) && (k2<(b+1)/2) &&
               (k3>=-(b/2)) && (k3<(b+1)/2)) /* left */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*3 + ((k2+b/2)*b+k3+b/2)*a + (k1+2*a);
	  else
	    return sum_N_B_less_r+N_B_r*3 + ((k1+2*a)*b+(k2+b/2))*b + (k3+b/2);
	}
      else if ((k2>=-2*a) && (k2<-a) && (k1>=-(b/2)) && (k1<(b+1)/2) &&
               (k3>=-(b/2)) && (k3<(b+1)/2)) /* front */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*4 + ((k1+b/2)*b+k3+b/2)*a + (k2+2*a);
	  else if (a==b)
	    return sum_N_B_less_r+N_B_r*4 + ((k1+b/2)*b+(k2+2*a))*a + (k3+b/2);
	  else
	    return sum_N_B_less_r+N_B_r*4 + ((k2+2*a)*b+(k1+b/2))*b + (k3+b/2);
	}
       else if ((k3>=-2*a) && (k3<-a) && (k1>=-(b/2)) && (k1<(b+1)/2) &&
                (k2>=-(b/2)) && (k2<(b+1)/2)) /* bottom */
	{
	  if(a>=b)
	    return sum_N_B_less_r+N_B_r*5 + ((k1+b/2)*b+k2+b/2)*a + (k3+2*a);
	  else
	    return sum_N_B_less_r+N_B_r*5 + ((k3+2*a)*b+(k1+b/2))*b + (k2+b/2);
	}

      sum_N_B_less_r+=6*N_B_r;
    } /* for(r) */

  return(-1);
}
#endif

#ifdef GAUSSIAN
static void init_index_sparse_to_full_3d(nsfft_plan *ths)
{
  int k1,k2,k3,k_s,r;
  int a,b;
  int N=X(exp2i)(ths->J+2);            /* length of the full grid           */
  int Nc=ths->center_nfft_plan->N[0];   /* length of the center block        */

  for (k_s=0, r=0; r<=(ths->J+1)/2; r++)
    {
      a=X(exp2i)(ths->J-r);
      b=X(exp2i)(r);

      /* right - rear - top - left - front - bottom */

      /* right */
      if(a>b)
	for(k2=-b/2;k2<(b+1)/2;k2++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k1=a; k1<2*a; k1++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k1=a; k1<2*a; k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;

      /* rear */
      if(a>b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k2=a; k2<2*a; k2++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else if(a==b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=a; k2<2*a; k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k2=a; k2<2*a; k2++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;

      /* top */
      if(a>=b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=a; k3<2*a; k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k3=a; k3<2*a; k3++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k2=-b/2;k2<(b+1)/2;k2++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;

      /* left */
      if(a>b)
	for(k2=-b/2;k2<(b+1)/2;k2++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k1=-2*a; k1<-a; k1++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k1=-2*a; k1<-a; k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;

      /* front */
      if(a>b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k2=-2*a; k2<-a; k2++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else if(a==b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=-2*a; k2<-a; k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k2=-2*a; k2<-a; k2++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;

      /* top */
      if(a>=b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=-2*a; k3<-a; k3++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k3=-2*a; k3<-a; k3++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k2=-b/2;k2<(b+1)/2;k2++,k_s++)
	      ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
    }

  /* center */
  for(k1=-Nc/2;k1<Nc/2;k1++)
    for(k2=-Nc/2;k2<Nc/2;k2++)
      for(k3=-Nc/2; k3<Nc/2; k3++,k_s++)
	ths->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
}
#endif

/* copies ths->f_hat to ths_plan->f_hat */
void nsfft_cp(nsfft_plan *ths, nfft_plan *ths_full_plan)
{
  int k;

  /* initialize f_hat with zero values */
  memset(ths_full_plan->f_hat, 0, ths_full_plan->N_total*sizeof(double _Complex));

   /* copy values at hyperbolic grid points */
  for(k=0;k<ths->N_total;k++)
    ths_full_plan->f_hat[ths->index_sparse_to_full[k]]=ths->f_hat[k];

  /* copy nodes */
  memcpy(ths_full_plan->x,ths->act_nfft_plan->x,ths->M_total*ths->d*sizeof(double));
}

#ifndef NSFTT_DISABLE_TEST
/* test copy_sparse_to_full */
static void test_copy_sparse_to_full_2d(nsfft_plan *ths, nfft_plan *ths_full_plan)
{
  int r;
  int k1, k2;
  int a,b;
  const int J=ths->J;   /* N=2^J                  */
  const int N=ths_full_plan->N[0];  /* size of full NFFT      */
  const int N_B=X(exp2i)(J);        /* size of small blocks   */

  /* copy sparse plan to full plan */
  nsfft_cp(ths, ths_full_plan);

  /* show blockwise f_hat */
  printf("f_hat blockwise\n");
  for (r=0; r<=(J+1)/2; r++){
    a=X(exp2i)(J-r); b=X(exp2i)(r);

    printf("top\n");
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        printf("(%1.1f,%1.1f) ", creal(ths->f_hat[(4*r+1)*N_B+ k1*b+k2]),
	                         cimag(ths->f_hat[(4*r+1)*N_B+ k1*b+k2]));
      }
      printf("\n");
    }

    printf("bottom\n");
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        printf("(%1.1f,%1.1f) ", creal(ths->f_hat[(4*r+3)*N_B+ k1*b+k2]),
                                 cimag(ths->f_hat[(4*r+3)*N_B+ k1*b+k2]));
      }
      printf("\n");
    }

    printf("right\n");
    for (k2=0; k2<b; k2++){
      for (k1=0; k1<a; k1++){
        printf("(%1.1f,%1.1f) ", creal(ths->f_hat[(4*r+0)*N_B+ k2*a+k1]),
                                 cimag(ths->f_hat[(4*r+0)*N_B+ k2*a+k1]));
      }
      printf("\n");
    }

    printf("left\n");
    for (k2=0; k2<b; k2++){
      for (k1=0; k1<a; k1++){
        printf("(%1.1f,%1.1f) ", creal(ths->f_hat[(4*r+2)*N_B+ k2*a+k1]),
	                         cimag(ths->f_hat[(4*r+2)*N_B+ k2*a+k1]));
      }
      printf("\n");
    }
  }

  return;
  /* show full f_hat */
  printf("full f_hat\n");
  for (k1=0;k1<N;k1++){
    for (k2=0;k2<N;k2++){
      printf("(%1.1f,%1.1f) ", creal(ths_full_plan->f_hat[k1*N+k2]),
                               cimag(ths_full_plan->f_hat[k1*N+k2]));
    }
    printf("\n");
  }
}
#endif

#ifndef NSFTT_DISABLE_TEST
static void test_sparse_to_full_2d(nsfft_plan* ths)
{
  int k_S,k1,k2;
  int N=X(exp2i)(ths->J+2);

  printf("N=%d\n\n",N);

  for(k1=0;k1<N;k1++)
    for(k2=0;k2<N;k2++)
      {
        k_S=index_full_to_sparse_2d(ths->J, k1*N+k2);
	if(k_S!=-1)
	  printf("(%+d, %+d)\t= %+d \t= %+d = %+d \n",k1-N/2,k2-N/2,
                 k1*N+k2, k_S, ths->index_sparse_to_full[k_S]);
      }
}
#endif

#ifndef NSFTT_DISABLE_TEST
static void test_sparse_to_full_3d(nsfft_plan* ths)
{
  int k_S,k1,k2,k3;
  int N=X(exp2i)(ths->J+2);

  printf("N=%d\n\n",N);

  for(k1=0;k1<N;k1++)
    for(k2=0;k2<N;k2++)
        for(k3=0;k3<N;k3++)
	  {
	    k_S=index_full_to_sparse_3d(ths->J, (k1*N+k2)*N+k3);
	    if(k_S!=-1)
	      printf("(%d, %d, %d)\t= %d \t= %d = %d \n",k1-N/2,k2-N/2,k3-N/2,
                     (k1*N+k2)*N+k3,k_S, ths->index_sparse_to_full[k_S]);
	  }
}
#endif


void nsfft_init_random_nodes_coeffs(nsfft_plan *ths)
{
  int j;

  /* init frequencies */
  nfft_vrand_unit_complex(ths->f_hat, ths->N_total);

  /* init nodes */
  nfft_vrand_shifted_unit_double(ths->act_nfft_plan->x, ths->d * ths->M_total);

  if(ths->d==2)
    for(j=0;j<ths->M_total;j++)
      {
        ths->x_transposed[2*j+0]=ths->act_nfft_plan->x[2*j+1];
        ths->x_transposed[2*j+1]=ths->act_nfft_plan->x[2*j+0];
      }
  else /* this->d==3 */
    for(j=0;j<ths->M_total;j++)
      {
        ths->x_102[3*j+0]=ths->act_nfft_plan->x[3*j+1];
        ths->x_102[3*j+1]=ths->act_nfft_plan->x[3*j+0];
        ths->x_102[3*j+2]=ths->act_nfft_plan->x[3*j+2];

        ths->x_201[3*j+0]=ths->act_nfft_plan->x[3*j+2];
        ths->x_201[3*j+1]=ths->act_nfft_plan->x[3*j+0];
        ths->x_201[3*j+2]=ths->act_nfft_plan->x[3*j+1];

        ths->x_120[3*j+0]=ths->act_nfft_plan->x[3*j+1];
        ths->x_120[3*j+1]=ths->act_nfft_plan->x[3*j+2];
        ths->x_120[3*j+2]=ths->act_nfft_plan->x[3*j+0];

        ths->x_021[3*j+0]=ths->act_nfft_plan->x[3*j+0];
        ths->x_021[3*j+1]=ths->act_nfft_plan->x[3*j+2];
        ths->x_021[3*j+2]=ths->act_nfft_plan->x[3*j+1];
      }
}

static void nsdft_trafo_2d(nsfft_plan *ths)
{
  int j,k_S,k_L,k0,k1;
  double omega;
  int N=X(exp2i)(ths->J+2);

  memset(ths->f,0,ths->M_total*sizeof(double _Complex));

  for(k_S=0;k_S<ths->N_total;k_S++)
    {
      k_L=ths->index_sparse_to_full[k_S];
      k0=k_L / N;
      k1=k_L % N;

      for(j=0;j<ths->M_total;j++)
	{
	  omega =
	    ((double)(k0 - N/2)) * ths->act_nfft_plan->x[2 * j + 0] +
	    ((double)(k1 - N/2)) * ths->act_nfft_plan->x[2 * j + 1];
          ths->f[j] += ths->f_hat[k_S] * cexp( - I*2*KPI*omega);
	}
    }
} /* void nsdft_trafo_2d */

static void nsdft_trafo_3d(nsfft_plan *ths)
{
  int j,k_S,k0,k1,k2;
  double omega;
  int N=X(exp2i)(ths->J+2);
  int k_L;

  memset(ths->f,0,ths->M_total*sizeof(double _Complex));

  for(k_S=0;k_S<ths->N_total;k_S++)
    {
      k_L=ths->index_sparse_to_full[k_S];

      k0=k_L/(N*N);
      k1=(k_L/N)%N;
      k2=k_L%N;

      for(j=0;j<ths->M_total;j++)
	{
	  omega =
	    ((double)(k0 - N/2)) * ths->act_nfft_plan->x[3 * j + 0] +
	    ((double)(k1 - N/2)) * ths->act_nfft_plan->x[3 * j + 1] +
	    ((double)(k2 - N/2)) * ths->act_nfft_plan->x[3 * j + 2];
          ths->f[j] += ths->f_hat[k_S] * cexp( - I*2*KPI*omega);
	}
    }
} /* void nsdft_trafo_3d */

void nsfft_trafo_direct(nsfft_plan *ths)
{
  if(ths->d==2)
    nsdft_trafo_2d(ths);
  else
    nsdft_trafo_3d(ths);
}

static void nsdft_adjoint_2d(nsfft_plan *ths)
{
  int j,k_S,k_L,k0,k1;
  double omega;
  int N=X(exp2i)(ths->J+2);

  memset(ths->f_hat,0,ths->N_total*sizeof(double _Complex));

  for(k_S=0;k_S<ths->N_total;k_S++)
    {
      k_L=ths->index_sparse_to_full[k_S];
      k0=k_L / N;
      k1=k_L % N;

      for(j=0;j<ths->M_total;j++)
	{
	  omega =
	    ((double)(k0 - N/2)) * ths->act_nfft_plan->x[2 * j + 0] +
	    ((double)(k1 - N/2)) * ths->act_nfft_plan->x[2 * j + 1];
          ths->f_hat[k_S] += ths->f[j] * cexp( + _Complex_I*2*KPI*omega);
	}
    }
} /* void nsdft_adjoint_2d */

static void nsdft_adjoint_3d(nsfft_plan *ths)
{
  int j,k_S,k0,k1,k2;
  double omega;
  int N=X(exp2i)(ths->J+2);
  int k_L;

  memset(ths->f_hat,0,ths->N_total*sizeof(double _Complex));

  for(k_S=0;k_S<ths->N_total;k_S++)
    {
      k_L=ths->index_sparse_to_full[k_S];

      k0=k_L/(N*N);
      k1=(k_L/N)%N;
      k2=k_L%N;

      for(j=0;j<ths->M_total;j++)
	{
	  omega =
	    ((double)(k0 - N/2)) * ths->act_nfft_plan->x[3 * j + 0] +
	    ((double)(k1 - N/2)) * ths->act_nfft_plan->x[3 * j + 1] +
	    ((double)(k2 - N/2)) * ths->act_nfft_plan->x[3 * j + 2];
          ths->f_hat[k_S] += ths->f[j] * cexp( + _Complex_I*2*KPI*omega);
	}
    }
} /* void nsdft_adjoint_3d */

void nsfft_adjoint_direct(nsfft_plan *ths)
{
  if(ths->d==2)
    nsdft_adjoint_2d(ths);
  else
    nsdft_adjoint_3d(ths);
}

static void nsfft_trafo_2d(nsfft_plan *ths)
{
  int r,rr,j;
  double temp;

  int M=ths->M_total;
  int J=ths->J;

  /* center */
  ths->center_nfft_plan->f_hat=ths->f_hat+4*((J+1)/2+1)*X(exp2i)(J);

  if (ths->center_nfft_plan->N[0]<=ths->center_nfft_plan->m)
    nfft_trafo_direct(ths->center_nfft_plan);
  else
    nfft_trafo(ths->center_nfft_plan);

  for (j=0; j<M; j++)
    ths->f[j] = ths->center_nfft_plan->f[j];

  for(rr=0;rr<=(J+1)/2;rr++)
    {
      r=MIN(rr,J-rr);
      ths->act_nfft_plan->my_fftw_plan1 = ths->set_fftw_plan1[r];
      ths->act_nfft_plan->N[0]=X(exp2i)(r); ths->act_nfft_plan->n[0]=ths->sigma*ths->act_nfft_plan->N[0];
      ths->act_nfft_plan->N[1]=X(exp2i)(J-r); ths->act_nfft_plan->n[1]=ths->sigma*ths->act_nfft_plan->N[1];

      /*printf("%d x %d\n",ths->act_nfft_plan->N[0],ths->act_nfft_plan->N[1]);*/

      temp=-3.0*KPI*X(exp2i)(J-rr);

      /* right */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+0)*X(exp2i)(J);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_trafo_direct(ths->act_nfft_plan);
	else
	  short_nfft_trafo_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      for (j=0; j<M; j++)
        ths->f[j] +=  ths->act_nfft_plan->f[j] *
                      cexp( + _Complex_I*temp*ths->act_nfft_plan->x[2*j+1]);

      /* top */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+1)*X(exp2i)(J);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_trafo_direct(ths->act_nfft_plan);
	else
	  short_nfft_trafo_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( + _Complex_I*temp*ths->act_nfft_plan->x[2*j+0]);

      /* left */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+2)*X(exp2i)(J);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_trafo_direct(ths->act_nfft_plan);
	else
	  short_nfft_trafo_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( - I*temp*ths->act_nfft_plan->x[2*j+1]);

      /* bottom */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+3)*X(exp2i)(J);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_trafo_direct(ths->act_nfft_plan);
	else
	  short_nfft_trafo_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( - I*temp*ths->act_nfft_plan->x[2*j+0]);
    } /* for(rr) */
} /* void nsfft_trafo_2d */

static void nsfft_adjoint_2d(nsfft_plan *ths)
{
  int r,rr,j;
  double temp;

  int M=ths->M_total;
  int J=ths->J;

  /* center */
  for (j=0; j<M; j++)
    ths->center_nfft_plan->f[j] = ths->f[j];

  ths->center_nfft_plan->f_hat=ths->f_hat+4*((J+1)/2+1)*X(exp2i)(J);

  if (ths->center_nfft_plan->N[0]<=ths->center_nfft_plan->m)
    nfft_adjoint_direct(ths->center_nfft_plan);
  else
    nfft_adjoint(ths->center_nfft_plan);

  for(rr=0;rr<=(J+1)/2;rr++)
    {
      r=MIN(rr,J-rr);
      ths->act_nfft_plan->my_fftw_plan2 = ths->set_fftw_plan2[r];
      ths->act_nfft_plan->N[0]=X(exp2i)(r); ths->act_nfft_plan->n[0]=ths->sigma*ths->act_nfft_plan->N[0];
      ths->act_nfft_plan->N[1]=X(exp2i)(J-r); ths->act_nfft_plan->n[1]=ths->sigma*ths->act_nfft_plan->N[1];

      /*printf("%d x %d\n",ths->act_nfft_plan->N[0],ths->act_nfft_plan->N[1]);*/

      temp=-3.0*KPI*X(exp2i)(J-rr);

      /* right */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+0)*X(exp2i)(J);

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( - I*temp*ths->act_nfft_plan->x[2*j+1]);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_adjoint_direct(ths->act_nfft_plan);
	else
	  short_nfft_adjoint_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      /* top */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+1)*X(exp2i)(J);

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( - I*temp*ths->act_nfft_plan->x[2*j+0]);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_adjoint_direct(ths->act_nfft_plan);
	else
	  short_nfft_adjoint_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      /* left */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+2)*X(exp2i)(J);

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( + _Complex_I*temp*ths->act_nfft_plan->x[2*j+1]);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_adjoint_direct(ths->act_nfft_plan);
	else
	  short_nfft_adjoint_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(r<rr)
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      /* bottom */
      ths->act_nfft_plan->f_hat=ths->f_hat+(4*rr+3)*X(exp2i)(J);

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( + _Complex_I*temp*ths->act_nfft_plan->x[2*j+0]);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  nfft_adjoint_direct(ths->act_nfft_plan);
	else
	  short_nfft_adjoint_2d(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if((r==rr)&&(J-rr!=rr))
	RSWAP(ths->act_nfft_plan->x,ths->x_transposed);
    } /* for(rr) */
} /* void nsfft_adjoint_2d */

static void nsfft_trafo_3d(nsfft_plan *ths)
{
  int r,rr,j;
  double temp;
  int sum_N_B_less_r,N_B_r,a,b;

  int M=ths->M_total;
  int J=ths->J;

  /* center */
  ths->center_nfft_plan->f_hat=ths->f_hat+6*X(exp2i)(J)*(X(exp2i)((J+1)/2+1)-1);

  if (ths->center_nfft_plan->N[0]<=ths->center_nfft_plan->m)
    nfft_trafo_direct(ths->center_nfft_plan);
  else
    nfft_trafo(ths->center_nfft_plan);

  for (j=0; j<M; j++)
    ths->f[j] = ths->center_nfft_plan->f[j];

  sum_N_B_less_r=0;
  for(rr=0;rr<=(J+1)/2;rr++)
    {
      a=X(exp2i)(J-rr);
      b=X(exp2i)(rr);

      N_B_r=a*b*b;

      r=MIN(rr,J-rr);
      ths->act_nfft_plan->my_fftw_plan1 = ths->set_fftw_plan1[rr];

      ths->act_nfft_plan->N[0]=X(exp2i)(r);
      if(a<b)
	ths->act_nfft_plan->N[1]=X(exp2i)(J-r);
      else
	ths->act_nfft_plan->N[1]=X(exp2i)(r);
      ths->act_nfft_plan->N[2]=X(exp2i)(J-r);

      /*printf("\n\n%d x %d x %d:\t",ths->act_nfft_plan->N[0],ths->act_nfft_plan->N[1],ths->act_nfft_plan->N[2]); fflush(stdout);*/

      ths->act_nfft_plan->N_total=ths->act_nfft_plan->N[0]*ths->act_nfft_plan->N[1]*ths->act_nfft_plan->N[2];
      ths->act_nfft_plan->n[0]=ths->sigma*ths->act_nfft_plan->N[0];
      ths->act_nfft_plan->n[1]=ths->sigma*ths->act_nfft_plan->N[1];
      ths->act_nfft_plan->n[2]=ths->sigma*ths->act_nfft_plan->N[2];
      ths->act_nfft_plan->n_total=ths->act_nfft_plan->n[0]*ths->act_nfft_plan->n[1]*ths->act_nfft_plan->n[2];

      /* only for right - rear - top */
      if((J==0)||((J==1)&&(rr==1)))
	temp=-2.0*KPI;
      else
	temp=-3.0*KPI*X(exp2i)(J-rr);

      /* right */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*0;

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_trafo_direct(ths->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( + _Complex_I*temp*ths->act_nfft_plan->x[3*j+0]);

      /* rear */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*1;

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_trafo_direct(ths->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( + _Complex_I*temp*ths->act_nfft_plan->x[3*j+1]);

      /* top */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*2;

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_trafo_direct(ths->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( + _Complex_I*temp*ths->act_nfft_plan->x[3*j+2]);

      /* only for left - front - bottom */
      if((J==0)||((J==1)&&(rr==1)))
	temp=-4.0*KPI;
      else
	temp=-3.0*KPI*X(exp2i)(J-rr);

      /* left */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*3;

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_trafo_direct(ths->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( - I*temp*ths->act_nfft_plan->x[3*j+0]);

      /* front */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*4;

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_trafo_direct(ths->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( - I*temp*ths->act_nfft_plan->x[3*j+1]);

      /* bottom */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*5;

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_trafo_direct(ths->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_trafo(ths->act_nfft_plan);

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      for (j=0; j<M; j++)
        ths->f[j] += ths->act_nfft_plan->f[j] *
                     cexp( - I*temp*ths->act_nfft_plan->x[3*j+2]);

      sum_N_B_less_r+=6*N_B_r;
    } /* for(rr) */
} /* void nsfft_trafo_3d */

static void nsfft_adjoint_3d(nsfft_plan *ths)
{
  int r,rr,j;
  double temp;
  int sum_N_B_less_r,N_B_r,a,b;

  int M=ths->M_total;
  int J=ths->J;

  /* center */
  for (j=0; j<M; j++)
    ths->center_nfft_plan->f[j] = ths->f[j];

  ths->center_nfft_plan->f_hat=ths->f_hat+6*X(exp2i)(J)*(X(exp2i)((J+1)/2+1)-1);

  if (ths->center_nfft_plan->N[0]<=ths->center_nfft_plan->m)
    nfft_adjoint_direct(ths->center_nfft_plan);
  else
    nfft_adjoint(ths->center_nfft_plan);

  sum_N_B_less_r=0;
  for(rr=0;rr<=(J+1)/2;rr++)
    {
      a=X(exp2i)(J-rr);
      b=X(exp2i)(rr);

      N_B_r=a*b*b;

      r=MIN(rr,J-rr);
      ths->act_nfft_plan->my_fftw_plan1 = ths->set_fftw_plan1[rr];
      ths->act_nfft_plan->my_fftw_plan2 = ths->set_fftw_plan2[rr];

      ths->act_nfft_plan->N[0]=X(exp2i)(r);
      if(a<b)
	ths->act_nfft_plan->N[1]=X(exp2i)(J-r);
      else
	ths->act_nfft_plan->N[1]=X(exp2i)(r);
      ths->act_nfft_plan->N[2]=X(exp2i)(J-r);

      /*printf("\n\n%d x %d x %d:\t",ths->act_nfft_plan->N[0],ths->act_nfft_plan->N[1],ths->act_nfft_plan->N[2]); fflush(stdout);*/

      ths->act_nfft_plan->N_total=ths->act_nfft_plan->N[0]*ths->act_nfft_plan->N[1]*ths->act_nfft_plan->N[2];
      ths->act_nfft_plan->n[0]=ths->sigma*ths->act_nfft_plan->N[0];
      ths->act_nfft_plan->n[1]=ths->sigma*ths->act_nfft_plan->N[1];
      ths->act_nfft_plan->n[2]=ths->sigma*ths->act_nfft_plan->N[2];
      ths->act_nfft_plan->n_total=ths->act_nfft_plan->n[0]*ths->act_nfft_plan->n[1]*ths->act_nfft_plan->n[2];

      /* only for right - rear - top */
      if((J==0)||((J==1)&&(rr==1)))
	temp=-2.0*KPI;
      else
	temp=-3.0*KPI*X(exp2i)(J-rr);

      /* right */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*0;

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( - I*temp*ths->act_nfft_plan->x[3*j+0]);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_adjoint_direct(ths->act_nfft_plan);
	  else
	    short_nfft_adjoint_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_adjoint_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      /* rear */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*1;

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( - I*temp*ths->act_nfft_plan->x[3*j+1]);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_adjoint_direct(ths->act_nfft_plan);
	  else
	    short_nfft_adjoint_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_adjoint_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      /* top */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*2;

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( - I*temp*ths->act_nfft_plan->x[3*j+2]);

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_adjoint_direct(ths->act_nfft_plan);
	  else
	    short_nfft_adjoint_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_adjoint_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      /* only for left - front - bottom */
      if((J==0)||((J==1)&&(rr==1)))
	temp=-4.0*KPI;
      else
	temp=-3.0*KPI*X(exp2i)(J-rr);

      /* left */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*3;

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( + _Complex_I*temp*ths->act_nfft_plan->x[3*j+0]);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_adjoint_direct(ths->act_nfft_plan);
	  else
	    short_nfft_adjoint_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_adjoint_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_120);

      /* front */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*4;

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( + _Complex_I*temp*ths->act_nfft_plan->x[3*j+1]);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_adjoint_direct(ths->act_nfft_plan);
	  else
	    short_nfft_adjoint_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_adjoint_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(a>b)
	RSWAP(ths->act_nfft_plan->x,ths->x_021);
      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_102);

      /* bottom */
      ths->act_nfft_plan->f_hat=ths->f_hat + sum_N_B_less_r + N_B_r*5;

      for (j=0; j<M; j++)
        ths->act_nfft_plan->f[j]= ths->f[j] *
                                  cexp( + _Complex_I*temp*ths->act_nfft_plan->x[3*j+2]);

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      if(ths->act_nfft_plan->N[0]<=ths->act_nfft_plan->m)
	if(ths->act_nfft_plan->N[1]<=ths->act_nfft_plan->m)
	  if(ths->act_nfft_plan->N[2]<=ths->act_nfft_plan->m)
	    nfft_adjoint_direct(ths->act_nfft_plan);
	  else
	    short_nfft_adjoint_3d_1(ths->act_nfft_plan,&(ths->set_nfft_plan_1d[r]));
	else
	  short_nfft_adjoint_3d_2(ths->act_nfft_plan,&(ths->set_nfft_plan_2d[r]));
      else
	nfft_adjoint(ths->act_nfft_plan);

      if(a<b)
	RSWAP(ths->act_nfft_plan->x,ths->x_201);

      sum_N_B_less_r+=6*N_B_r;
    } /* for(rr) */
} /* void nsfft_adjoint_3d */

void nsfft_trafo(nsfft_plan *ths)
{
  if(ths->d==2)
    nsfft_trafo_2d(ths);
  else
    nsfft_trafo_3d(ths);
}

void nsfft_adjoint(nsfft_plan *ths)
{
  if(ths->d==2)
    nsfft_adjoint_2d(ths);
  else
    nsfft_adjoint_3d(ths);
}


/*========================================================*/
/* J >1, no precomputation at all!! */
#ifdef GAUSSIAN
static void nsfft_init_2d(nsfft_plan *ths, int J, int M, int m, unsigned snfft_flags)
{
  int r;
  int N[2];
  int n[2];

  ths->flags=snfft_flags;
  ths->sigma=2;
  ths->J=J;
  ths->M_total=M;
  ths->N_total=(J+4)*X(exp2i)(J+1);

  /* memory allocation */
  ths->f = (double _Complex *)nfft_malloc(M*sizeof(double _Complex));
  ths->f_hat = (double _Complex *)nfft_malloc(ths->N_total*sizeof(double _Complex));
  ths->x_transposed= (double*)nfft_malloc(2*M*sizeof(double));

  ths->act_nfft_plan = (nfft_plan*)nfft_malloc(sizeof(nfft_plan));
  ths->center_nfft_plan = (nfft_plan*)nfft_malloc(sizeof(nfft_plan));

  ths->set_fftw_plan1=(fftw_plan*) nfft_malloc((J/2+1)*sizeof(fftw_plan));
  ths->set_fftw_plan2=(fftw_plan*) nfft_malloc((J/2+1)*sizeof(fftw_plan));

  ths->set_nfft_plan_1d = (nfft_plan*) nfft_malloc((X(log2i)(m)+1)*(sizeof(nfft_plan)));

  /* planning the small nffts */
  /* r=0 */
  N[0]=1;            n[0]=ths->sigma*N[0];
  N[1]=X(exp2i)(J); n[1]=ths->sigma*N[1];

  nfft_init_guru(ths->act_nfft_plan,2,N,M,n,m, FG_PSI| MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);

  if(ths->act_nfft_plan->flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(ths->act_nfft_plan);

  ths->set_fftw_plan1[0]=ths->act_nfft_plan->my_fftw_plan1;
  ths->set_fftw_plan2[0]=ths->act_nfft_plan->my_fftw_plan2;

  for(r=1;r<=J/2;r++)
    {
      N[0]=X(exp2i)(r);   n[0]=ths->sigma*N[0];
      N[1]=X(exp2i)(J-r); n[1]=ths->sigma*N[1];
      ths->set_fftw_plan1[r] =
	fftw_plan_dft(2, n, ths->act_nfft_plan->g1, ths->act_nfft_plan->g2,
		      FFTW_FORWARD, ths->act_nfft_plan->fftw_flags);

      ths->set_fftw_plan2[r] =
	fftw_plan_dft(2, n, ths->act_nfft_plan->g2, ths->act_nfft_plan->g1,
		      FFTW_BACKWARD, ths->act_nfft_plan->fftw_flags);
    }

  /* planning the 1d nffts */
  for(r=0;r<=X(log2i)(m);r++)
    {
      N[0]=X(exp2i)(J-r); n[0]=ths->sigma*N[0]; /* ==N[1] of the 2 dimensional plan */

      nfft_init_guru(&(ths->set_nfft_plan_1d[r]),1,N,M,n,m, MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
      ths->set_nfft_plan_1d[r].flags = ths->set_nfft_plan_1d[r].flags | FG_PSI;
      ths->set_nfft_plan_1d[r].K=ths->act_nfft_plan->K;
      ths->set_nfft_plan_1d[r].psi=ths->act_nfft_plan->psi;
    }

  /* center plan */
  /* J/2 == floor(((double)J) / 2.0) */
  N[0]=X(exp2i)(J/2+1); n[0]=ths->sigma*N[0];
  N[1]=X(exp2i)(J/2+1); n[1]=ths->sigma*N[1];
  nfft_init_guru(ths->center_nfft_plan,2,N,M,n, m, MALLOC_F| FFTW_INIT,
                     FFTW_MEASURE);
  ths->center_nfft_plan->x= ths->act_nfft_plan->x;
  ths->center_nfft_plan->flags = ths->center_nfft_plan->flags|
                                      FG_PSI;
  ths->center_nfft_plan->K=ths->act_nfft_plan->K;
  ths->center_nfft_plan->psi=ths->act_nfft_plan->psi;

  if(ths->flags & NSDFT)
    {
      ths->index_sparse_to_full=(int*)nfft_malloc(ths->N_total*sizeof(int));
      init_index_sparse_to_full_2d(ths);
    }
}
#endif

/*========================================================*/
/* J >1, no precomputation at all!! */
#ifdef GAUSSIAN
static void nsfft_init_3d(nsfft_plan *ths, int J, int M, int m, unsigned snfft_flags)
{
  int r,rr,a,b;
  int N[3];
  int n[3];

  ths->flags=snfft_flags;
  ths->sigma=2;
  ths->J=J;
  ths->M_total=M;
  ths->N_total=6*X(exp2i)(J)*(X(exp2i)((J+1)/2+1)-1)+X(exp2i)(3*(J/2+1));

  /* memory allocation */
  ths->f =     (double _Complex *)nfft_malloc(M*sizeof(double _Complex));
  ths->f_hat = (double _Complex *)nfft_malloc(ths->N_total*sizeof(double _Complex));

  ths->x_102= (double*)nfft_malloc(3*M*sizeof(double));
  ths->x_201= (double*)nfft_malloc(3*M*sizeof(double));
  ths->x_120= (double*)nfft_malloc(3*M*sizeof(double));
  ths->x_021= (double*)nfft_malloc(3*M*sizeof(double));

  ths->act_nfft_plan = (nfft_plan*)nfft_malloc(sizeof(nfft_plan));
  ths->center_nfft_plan = (nfft_plan*)nfft_malloc(sizeof(nfft_plan));

  ths->set_fftw_plan1=(fftw_plan*) nfft_malloc(((J+1)/2+1)*sizeof(fftw_plan));
  ths->set_fftw_plan2=(fftw_plan*) nfft_malloc(((J+1)/2+1)*sizeof(fftw_plan));

  ths->set_nfft_plan_1d = (nfft_plan*) nfft_malloc((X(log2i)(m)+1)*(sizeof(nfft_plan)));
  ths->set_nfft_plan_2d = (nfft_plan*) nfft_malloc((X(log2i)(m)+1)*(sizeof(nfft_plan)));

  /* planning the small nffts */
  /* r=0 */
  N[0]=1;            n[0]=ths->sigma*N[0];
  N[1]=1;            n[1]=ths->sigma*N[1];
  N[2]=X(exp2i)(J); n[2]=ths->sigma*N[2];

  nfft_init_guru(ths->act_nfft_plan,3,N,M,n,m, FG_PSI| MALLOC_X| MALLOC_F, FFTW_MEASURE);

  if(ths->act_nfft_plan->flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(ths->act_nfft_plan);

  /* malloc g1, g2 for maximal size */
  ths->act_nfft_plan->g1 = nfft_malloc(ths->sigma*ths->sigma*ths->sigma*X(exp2i)(J+(J+1)/2)*sizeof(double _Complex));
  ths->act_nfft_plan->g2 = nfft_malloc(ths->sigma*ths->sigma*ths->sigma*X(exp2i)(J+(J+1)/2)*sizeof(double _Complex));

  ths->act_nfft_plan->my_fftw_plan1 =
    fftw_plan_dft(3, n, ths->act_nfft_plan->g1, ths->act_nfft_plan->g2,
		  FFTW_FORWARD, ths->act_nfft_plan->fftw_flags);
  ths->act_nfft_plan->my_fftw_plan2 =
    fftw_plan_dft(3, n, ths->act_nfft_plan->g2, ths->act_nfft_plan->g1,
		  FFTW_BACKWARD, ths->act_nfft_plan->fftw_flags);

  ths->set_fftw_plan1[0]=ths->act_nfft_plan->my_fftw_plan1;
  ths->set_fftw_plan2[0]=ths->act_nfft_plan->my_fftw_plan2;

  for(rr=1;rr<=(J+1)/2;rr++)
    {
      a=X(exp2i)(J-rr);
      b=X(exp2i)(rr);

      r=MIN(rr,J-rr);

      n[0]=ths->sigma*X(exp2i)(r);
      if(a<b)
	n[1]=ths->sigma*X(exp2i)(J-r);
      else
	n[1]=ths->sigma*X(exp2i)(r);
      n[2]=ths->sigma*X(exp2i)(J-r);

      ths->set_fftw_plan1[rr] =
	fftw_plan_dft(3, n, ths->act_nfft_plan->g1, ths->act_nfft_plan->g2,
		      FFTW_FORWARD, ths->act_nfft_plan->fftw_flags);
      ths->set_fftw_plan2[rr] =
	fftw_plan_dft(3, n, ths->act_nfft_plan->g2, ths->act_nfft_plan->g1,
		      FFTW_BACKWARD, ths->act_nfft_plan->fftw_flags);
    }

  /* planning the 1d nffts */
  for(r=0;r<=X(log2i)(m);r++)
    {
      N[0]=X(exp2i)(J-r); n[0]=ths->sigma*N[0];
      N[1]=X(exp2i)(J-r); n[1]=ths->sigma*N[1];

      if(N[0]>m)
	{
	  nfft_init_guru(&(ths->set_nfft_plan_1d[r]),1,N,M,n,m, MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
	  ths->set_nfft_plan_1d[r].flags = ths->set_nfft_plan_1d[r].flags | FG_PSI;
	  ths->set_nfft_plan_1d[r].K=ths->act_nfft_plan->K;
	  ths->set_nfft_plan_1d[r].psi=ths->act_nfft_plan->psi;
	  nfft_init_guru(&(ths->set_nfft_plan_2d[r]),2,N,M,n,m, MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
	  ths->set_nfft_plan_2d[r].flags = ths->set_nfft_plan_2d[r].flags | FG_PSI;
	  ths->set_nfft_plan_2d[r].K=ths->act_nfft_plan->K;
	  ths->set_nfft_plan_2d[r].psi=ths->act_nfft_plan->psi;
	}
    }

  /* center plan */
  /* J/2 == floor(((double)J) / 2.0) */
  N[0]=X(exp2i)(J/2+1); n[0]=ths->sigma*N[0];
  N[1]=X(exp2i)(J/2+1); n[1]=ths->sigma*N[1];
  N[2]=X(exp2i)(J/2+1); n[2]=ths->sigma*N[2];
  nfft_init_guru(ths->center_nfft_plan,3,N,M,n, m, MALLOC_F| FFTW_INIT,
                     FFTW_MEASURE);
  ths->center_nfft_plan->x= ths->act_nfft_plan->x;
  ths->center_nfft_plan->flags = ths->center_nfft_plan->flags|
                                      FG_PSI;
  ths->center_nfft_plan->K=ths->act_nfft_plan->K;
  ths->center_nfft_plan->psi=ths->act_nfft_plan->psi;

  if(ths->flags & NSDFT)
    {
      ths->index_sparse_to_full=(int*)nfft_malloc(ths->N_total*sizeof(int));
      init_index_sparse_to_full_3d(ths);
    }
}
#endif

#ifdef GAUSSIAN
void nsfft_init(nsfft_plan *ths, int d, int J, int M, int m, unsigned flags)
{
  ths->d=d;

  if(ths->d==2)
    nsfft_init_2d(ths, J, M, m, flags);
  else
    nsfft_init_3d(ths, J, M, m, flags);


  ths->mv_trafo = (void (*) (void* ))nsfft_trafo;
  ths->mv_adjoint = (void (*) (void* ))nsfft_adjoint;
}
#else
void nsfft_init(nsfft_plan *ths, int d, int J, int M, int m, unsigned flags)
{
  UNUSED(ths);
  UNUSED(d);
  UNUSED(J);
  UNUSED(M);
  UNUSED(m);
  UNUSED(flags);
  fprintf(stderr,
	  "\nError in kernel/nsfft_init: require GAUSSIAN window function\n");
}
#endif

static void nsfft_finalize_2d(nsfft_plan *ths)
{
  int r;

  if(ths->flags & NSDFT)
    nfft_free(ths->index_sparse_to_full);

  /* center plan */
  ths->center_nfft_plan->flags = ths->center_nfft_plan->flags ^ FG_PSI;
  nfft_finalize(ths->center_nfft_plan);

  /* the 1d nffts */
  for(r=0;r<=X(log2i)(ths->act_nfft_plan->m);r++)
    {
      ths->set_nfft_plan_1d[r].flags =
        ths->set_nfft_plan_1d[r].flags ^ FG_PSI;
      nfft_finalize(&(ths->set_nfft_plan_1d[r]));
    }

  /* finalize the small nffts */
  ths->act_nfft_plan->my_fftw_plan2=ths->set_fftw_plan2[0];
  ths->act_nfft_plan->my_fftw_plan1=ths->set_fftw_plan1[0];

  for(r=1;r<=ths->J/2;r++)
    {
      fftw_destroy_plan(ths->set_fftw_plan2[r]);
      fftw_destroy_plan(ths->set_fftw_plan1[r]);
    }

  /* r=0 */
  nfft_finalize(ths->act_nfft_plan);

  nfft_free(ths->set_nfft_plan_1d);

  nfft_free(ths->set_fftw_plan2);
  nfft_free(ths->set_fftw_plan1);

  nfft_free(ths->x_transposed);

  nfft_free(ths->f_hat);
  nfft_free(ths->f);
}

static void nsfft_finalize_3d(nsfft_plan *ths)
{
  int r;

  if(ths->flags & NSDFT)
    nfft_free(ths->index_sparse_to_full);

  /* center plan */
  ths->center_nfft_plan->flags = ths->center_nfft_plan->flags ^ FG_PSI;
  nfft_finalize(ths->center_nfft_plan);

  /* the 1d and 2d nffts */
  for(r=0;r<=X(log2i)(ths->act_nfft_plan->m);r++)
    {
      if(X(exp2i)(ths->J-r)>ths->act_nfft_plan->m)
	{
	  ths->set_nfft_plan_2d[r].flags = ths->set_nfft_plan_2d[r].flags ^ FG_PSI;
	  nfft_finalize(&(ths->set_nfft_plan_2d[r]));
	  ths->set_nfft_plan_1d[r].flags = ths->set_nfft_plan_1d[r].flags ^ FG_PSI;
	  nfft_finalize(&(ths->set_nfft_plan_1d[r]));
	}
    }

  /* finalize the small nffts */
  ths->act_nfft_plan->my_fftw_plan2=ths->set_fftw_plan2[0];
  ths->act_nfft_plan->my_fftw_plan1=ths->set_fftw_plan1[0];

  for(r=1;r<=(ths->J+1)/2;r++)
    {
      fftw_destroy_plan(ths->set_fftw_plan2[r]);
      fftw_destroy_plan(ths->set_fftw_plan1[r]);
    }

  /* r=0 */
  nfft_finalize(ths->act_nfft_plan);

  nfft_free(ths->set_nfft_plan_1d);
  nfft_free(ths->set_nfft_plan_2d);

  nfft_free(ths->set_fftw_plan2);
  nfft_free(ths->set_fftw_plan1);

  nfft_free(ths->x_102);
  nfft_free(ths->x_201);
  nfft_free(ths->x_120);
  nfft_free(ths->x_021);

  nfft_free(ths->f_hat);
  nfft_free(ths->f);
}

void nsfft_finalize(nsfft_plan *ths)
{
  if(ths->d==2)
    nsfft_finalize_2d(ths);
  else
    nsfft_finalize_3d(ths);
}
