#include "sparse_nfft.h"

int total_used_memory()
{
  struct mallinfo m;
  m=mallinfo();
  return m.hblkhd + m.uordblks;
}

int int_2_pow(int a)
{
  return (1U<< a);
}

int ld(int m)
{
  int l=0;
  int mm=m;
  
  while(mm>0)
    {
      mm=(mm>>1);
      l++;
    }
  return (l-1);
}

int index_sparse_to_full_direct(int J, int k)
{
    int N=int_2_pow(J+2);               /* number of full coeffs            */
    int N_B=int_2_pow(J);               /* number in each sparse block      */

    int j=k/N_B;                        /* consecutive number of Block      */
    int r=j/4;                          /* level of block                   */

    int i, o, a, b,s,l,m1,m2;
    int k1,k2;

    if (k>=(J+4)*int_2_pow(J+1))
      {
	printf("Fehler!\n");
	return(-1);
      }
    else
      {
	if (r>(J+1)/2)                      /* center block                     */
	  {
	    i=k-4*((J+1)/2+1)*N_B;
	    a=int_2_pow(J/2+1);
	    m1=i/a;
	    m2=i%a;
	    k1=N/2-a/2+m1;
	    k2=N/2-a/2+m2;
	  }
	else                                /* no center block                  */
	  {
	    i=k-j*N_B;                      /* index in specific block          */
	    o=j%4;                          /* kind of specific block           */
	    a=int_2_pow(r);
	    b=int_2_pow(J-r);
	    l=MAX(a,b);                     /* long dimension of block     */
	    s=MIN(a,b);                     /* short dimension of block    */
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

inline int index_sparse_to_full(sparse_nfft_plan* this, int k)
{
  /* only by lookup table */
  if( k < this->N_S)
    return this->index_sparse_to_full[k];
  else
    return -1;
}
 
int index_full_to_sparse(int J, int k)
{
    int N=int_2_pow(J+2);  /* number of full coeffs       */
    int N_B=int_2_pow(J);  /* number in each sparse block */

    int k1=k/N-N/2;        /* coordinates in the full grid */
    int k2=k%N-N/2;        /* k1: row, k2: column          */

    int r,a,b;

    a=int_2_pow(J/2+1);

    if ( (k1>=-(a/2)) && (k1<a/2) && (k2>=(-a/2)) && (k2<a/2) )
      {
	return(4*((J+1)/2+1)*N_B+(k1+a/2)*a+(k2+a/2));
      }

    for (r=0; r<=(J+1)/2; r++)
      {
	b=int_2_pow(r);
	a=int_2_pow(J-r);
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

void init_index_sparse_to_full(sparse_nfft_plan* this)
{
  int k_S;

  for (k_S=0; k_S<this->N_S; k_S++)
    this->index_sparse_to_full[k_S]=index_sparse_to_full_direct(this->J, k_S);
}

inline int index_sparse_to_full_3d(sparse_nfft_plan_3d* this, int k)
{
  /* only by lookup table */
  if( k < this->N_S)
    return this->index_sparse_to_full[k];
  else
    return -1;
}

int index_full_to_sparse_3d(int J, int k)
{
  int N=int_2_pow(J+2);                 /* length of the full grid            */
  int N_B_r;                            /* size of a sparse block in level r  */
  int sum_N_B_less_r;                       /* sum N_B_r                          */

  int r,a,b;

  int k3=(k%N)-N/2;                       /* coordinates in the full grid       */
  int k2=((k/N)%N)-N/2;
  int k1=k/(N*N)-N/2;
    
  a=int_2_pow(J/2+1);                   /* length of center block             */

  if((k1>=-(a/2)) && (k1<a/2) && (k2>=(-a/2)) && (k2<a/2) && (k3>=(-a/2)) && (k3<a/2))  
    {
      return(6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1)+((k1+a/2)*a+(k2+a/2))*a+(k3+a/2));
    }

  sum_N_B_less_r=0;
  for (r=0; r<=(J+1)/2; r++)
    {
      a=int_2_pow(J-r);
      b=int_2_pow(r);

      N_B_r=a*b*b;

      /* right - rear - top - left - front - bottom */
      if ((k1>=a) && (k1<2*a) && (k2>=-(b/2)) && (k2<(b+1)/2) && (k3>=-(b/2)) && (k3<(b+1)/2)) /* right */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*0 + ((k2+b/2)*b+k3+b/2)*a + (k1-a);
	  else
	    return sum_N_B_less_r+N_B_r*0 + ((k1-a)*b+(k2+b/2))*b + (k3+b/2);
	}
      else if ((k2>=a) && (k2<2*a) && (k1>=-(b/2)) && (k1<(b+1)/2) && (k3>=-(b/2)) && (k3<(b+1)/2)) /* rear */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*1 + ((k1+b/2)*b+k3+b/2)*a + (k2-a);
	  else if (a==b)
	    return sum_N_B_less_r+N_B_r*1 + ((k1+b/2)*b+(k2-a))*a + (k3+b/2);
	  else
	    return sum_N_B_less_r+N_B_r*1 + ((k2-a)*b+(k1+b/2))*b + (k3+b/2);
	}
       else if ((k3>=a) && (k3<2*a) && (k1>=-(b/2)) && (k1<(b+1)/2) && (k2>=-(b/2)) && (k2<(b+1)/2)) /* top */
	{
	  if(a>=b)
	    return sum_N_B_less_r+N_B_r*2 + ((k1+b/2)*b+k2+b/2)*a + (k3-a);
	  else
	    return sum_N_B_less_r+N_B_r*2 + ((k3-a)*b+(k1+b/2))*b + (k2+b/2);
	}

      else if ((k1>=-2*a) && (k1<-a) && (k2>=-(b/2)) && (k2<(b+1)/2) && (k3>=-(b/2)) && (k3<(b+1)/2)) /* left */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*3 + ((k2+b/2)*b+k3+b/2)*a + (k1+2*a);
	  else
	    return sum_N_B_less_r+N_B_r*3 + ((k1+2*a)*b+(k2+b/2))*b + (k3+b/2);
	}
      else if ((k2>=-2*a) && (k2<-a) && (k1>=-(b/2)) && (k1<(b+1)/2) && (k3>=-(b/2)) && (k3<(b+1)/2)) /* front */
	{
	  if(a>b)
	    return sum_N_B_less_r+N_B_r*4 + ((k1+b/2)*b+k3+b/2)*a + (k2+2*a);
	  else if (a==b)
	    return sum_N_B_less_r+N_B_r*4 + ((k1+b/2)*b+(k2+2*a))*a + (k3+b/2);
	  else
	    return sum_N_B_less_r+N_B_r*4 + ((k2+2*a)*b+(k1+b/2))*b + (k3+b/2);
	}
       else if ((k3>=-2*a) && (k3<-a) && (k1>=-(b/2)) && (k1<(b+1)/2) && (k2>=-(b/2)) && (k2<(b+1)/2)) /* bottom */
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

void init_index_sparse_to_full_3d(sparse_nfft_plan_3d* this)
{
  int k1,k2,k3,k_s,r;
  int a,b;
  int N=int_2_pow(this->J+2);           /* length of the full grid            */
  int Nc=this->center_nfft_plan->N[0];  /* length of the center block         */

  for (k_s=0, r=0; r<=(this->J+1)/2; r++)
    {
      a=int_2_pow(this->J-r);
      b=int_2_pow(r);

      /* right - rear - top - left - front - bottom */

      /* right */
      if(a>b)
	for(k2=-b/2;k2<(b+1)/2;k2++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k1=a; k1<2*a; k1++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k1=a; k1<2*a; k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
	
      /* rear */
      if(a>b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k2=a; k2<2*a; k2++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else if(a==b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=a; k2<2*a; k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k2=a; k2<2*a; k2++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      
      /* top */
      if(a>=b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=a; k3<2*a; k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k3=a; k3<2*a; k3++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k2=-b/2;k2<(b+1)/2;k2++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
	
      /* left */
      if(a>b)
	for(k2=-b/2;k2<(b+1)/2;k2++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k1=-2*a; k1<-a; k1++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k1=-2*a; k1<-a; k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
	
      /* front */
      if(a>b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k3=-b/2;k3<(b+1)/2;k3++)
	    for(k2=-2*a; k2<-a; k2++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else if(a==b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=-2*a; k2<-a; k2++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k2=-2*a; k2<-a; k2++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k3=-b/2;k3<(b+1)/2;k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      
      /* top */
      if(a>=b)
	for(k1=-b/2;k1<(b+1)/2;k1++)
	  for(k2=-b/2;k2<(b+1)/2;k2++)
	    for(k3=-2*a; k3<-a; k3++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
      else
	for(k3=-2*a; k3<-a; k3++)
	  for(k1=-b/2;k1<(b+1)/2;k1++)
	    for(k2=-b/2;k2<(b+1)/2;k2++,k_s++)
	      this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
    }

  /* center */
  for(k1=-Nc/2;k1<Nc/2;k1++)
    for(k2=-Nc/2;k2<Nc/2;k2++)
      for(k3=-Nc/2; k3<Nc/2; k3++,k_s++)
	this->index_sparse_to_full[k_s]=((k1+N/2)*N+k2+N/2)*N+k3+N/2;
}

/* copies this->f_hat to this_plan->f_hat */
void copy_sparse_to_full_3d(sparse_nfft_plan_3d *this, nfft_plan *this_full_plan)
{
  int k_S,k; 
  int N=this_full_plan->N[0];     /* length of the full grid            */
  
  /* initialize f_hat with zero values */
  memset(this_full_plan->f_hat, 0, N*N*N*sizeof(fftw_complex));
  
   /* copy values at hyperbolic grid points */
  for(k_S=0;k_S<this->N_S;k_S++)
    {
      k=(int)this->index_sparse_to_full[k_S];
      this_full_plan->f_hat[k][0]=this->f_hat[k_S][0];
      this_full_plan->f_hat[k][1]=this->f_hat[k_S][1];
    }
  
  /* copy nodes */
  memcpy(this_full_plan->x,this->act_nfft_plan->x,this->M*3*sizeof(double));
}

void test_sparse_to_full_3d(sparse_nfft_plan_3d* this)
{
  int k_S,k1,k2,k3;
  int N=int_2_pow(this->J+2);

  printf("N=%d\n\n",N);

  for(k1=0;k1<N;k1++)
    for(k2=0;k2<N;k2++)
        for(k3=0;k3<N;k3++)
	  {
	    k_S=index_full_to_sparse_3d(this->J, (k1*N+k2)*N+k3);
	    if(k_S!=-1)
	      printf("(%d, %d, %d)\t= %d \t= %d = %d \n",k1-N/2,k2-N/2,k3-N/2,(k1*N+k2)*N+k3,k_S, this->index_sparse_to_full[k_S]);
	  }
}

/* copies this->f_hat to this_plan->f_hat */
void copy_sparse_to_full(sparse_nfft_plan *this, nfft_plan *this_full_plan)
{
  int k_S,k; 
  const int N=this_full_plan->N[0];  /* size of full NFFT      */
  

  /* initialize f_hat with zero values */
  memset(this_full_plan->f_hat, 0, N*N*sizeof(fftw_complex));
  
   /* copy values at hyperbolic grid points */
  for(k_S=0;k_S<this->N_S;k_S++)
    {
      k=(int)this->index_sparse_to_full[k_S];
      this_full_plan->f_hat[k][0]=this->f_hat[k_S][0];
      this_full_plan->f_hat[k][1]=this->f_hat[k_S][1];
    }
  
  /* copy nodes */
  memcpy(this_full_plan->x,this->act_nfft_plan->x,this->M*2*sizeof(double));
}

/* test copy_sparse_to_full */
void test_copy_sparse_to_full(sparse_nfft_plan *this, nfft_plan *this_full_plan)
{
  int r;
  int k1, k2;
  int a,b;
  const int J=this->J;   /* N=2^J                  */
  const int N=this_full_plan->N[0];  /* size of full NFFT      */
  const int N_B=int_2_pow(J);        /* size of small blocks   */

  /* copy sparse plan to full plan */
  copy_sparse_to_full(this, this_full_plan);

  /* show blockwise f_hat */
  printf("f_hat blockwise\n");
  for (r=0; r<=(J+1)/2; r++){
    a=int_2_pow(J-r); b=int_2_pow(r);

    printf("top\n");
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        printf("(%1.1f,%1.1f) ", this->f_hat[(4*r+1)*N_B+ k1*b+k2 ][0]
	                       , this->f_hat[(4*r+1)*N_B+ k1*b+k2 ][1]);
      }
      printf("\n");
    }

    printf("bottom\n");
    for (k1=0; k1<a; k1++){
      for (k2=0; k2<b; k2++){
        printf("(%1.1f,%1.1f) ", this->f_hat[(4*r+3)*N_B+ k1*b+k2 ][0]
                               , this->f_hat[(4*r+3)*N_B+ k1*b+k2 ][1]);
      }
      printf("\n");
    }

    printf("right\n");
    for (k2=0; k2<b; k2++){
      for (k1=0; k1<a; k1++){
        printf("(%1.1f,%1.1f) ", this->f_hat[(4*r+0)*N_B+ k2*a+k1 ][0]
                               , this->f_hat[(4*r+0)*N_B+ k2*a+k1 ][1]);
      }
      printf("\n");
    }

    printf("left\n");
    for (k2=0; k2<b; k2++){
      for (k1=0; k1<a; k1++){
        printf("(%1.1f,%1.1f) ", this->f_hat[(4*r+2)*N_B+ k2*a+k1 ][0]
	                       , this->f_hat[(4*r+2)*N_B+ k2*a+k1 ][1]);
      }
      printf("\n");
    }
  }

  return;
  /* show full f_hat */
  printf("full f_hat\n");
  for (k1=0;k1<N;k1++){
    for (k2=0;k2<N;k2++){
      printf("(%1.1f,%1.1f) ", this_full_plan->f_hat[k1*N+k2][0],this_full_plan->f_hat[k1*N+k2][1]);
    }
    printf("\n");
  }
}

void sparse_init_random_nodes_coeffs(sparse_nfft_plan *this)
{
  int k_S,j;

  /* init frequencies */
  for(k_S=0;k_S<this->N_S;k_S++)
    {
      this->f_hat[k_S][0]=(double)rand()/RAND_MAX;
      this->f_hat[k_S][1]=(double)rand()/RAND_MAX;
    }
 
  /* init nodes */
  for(j=0;j<this->M;j++) 
    {
      this->act_nfft_plan->x[2*j+0]=((double)rand())/RAND_MAX-0.5;
      this->act_nfft_plan->x[2*j+1]=((double)rand())/RAND_MAX-0.5;
      this->x_transposed[2*j+0]=this->act_nfft_plan->x[2*j+1];
      this->x_transposed[2*j+1]=this->act_nfft_plan->x[2*j+0];
    }
}

void sparse_init_random_nodes_coeffs_3d(sparse_nfft_plan_3d *this)
{
  int k_S,j;

  /* init frequencies */
  for(k_S=0;k_S<this->N_S;k_S++)
    {
      this->f_hat[k_S][0]=(double)rand()/RAND_MAX;
      this->f_hat[k_S][1]=(double)rand()/RAND_MAX;
    }

  //this->f_hat[672][0]=1;
 
  /* init nodes */
  for(j=0;j<this->M;j++) 
    {
      this->act_nfft_plan->x[3*j+0]=((double)rand())/RAND_MAX-0.5;
      this->act_nfft_plan->x[3*j+1]=((double)rand())/RAND_MAX-0.5;
      this->act_nfft_plan->x[3*j+2]=((double)rand())/RAND_MAX-0.5;
     
      this->x_102[3*j+0]=this->act_nfft_plan->x[3*j+1];
      this->x_102[3*j+1]=this->act_nfft_plan->x[3*j+0];	
      this->x_102[3*j+2]=this->act_nfft_plan->x[3*j+2];

      this->x_201[3*j+0]=this->act_nfft_plan->x[3*j+2];
      this->x_201[3*j+1]=this->act_nfft_plan->x[3*j+0];	
      this->x_201[3*j+2]=this->act_nfft_plan->x[3*j+1];

      this->x_120[3*j+0]=this->act_nfft_plan->x[3*j+1];
      this->x_120[3*j+1]=this->act_nfft_plan->x[3*j+2];	
      this->x_120[3*j+2]=this->act_nfft_plan->x[3*j+0];

      this->x_021[3*j+0]=this->act_nfft_plan->x[3*j+0];
      this->x_021[3*j+1]=this->act_nfft_plan->x[3*j+2];	
      this->x_021[3*j+2]=this->act_nfft_plan->x[3*j+1];
    }
}

void sparse_ndft_trafo(sparse_nfft_plan *this)
{
  int j,k_S,k_L,k0,k1;
  double omega,cosx,sinx;
  int N=int_2_pow(this->J+2);

  memset(this->f,0,this->M*sizeof(fftw_complex));

  for(k_S=0;k_S<this->N_S;k_S++)
    {
      k_L=this->index_sparse_to_full[k_S];
      k0=k_L / N;
      k1=k_L % N;
      
      for(j=0;j<this->M;j++)
	{
	  omega =
	    ((double)(k0 - N/2)) * this->act_nfft_plan->x[2 * j + 0] + 
	    ((double)(k1 - N/2)) * this->act_nfft_plan->x[2 * j + 1];
	  cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);

	  this->f[j][0] += (this->f_hat[k_S][0]*cosx + this->f_hat[k_S][1]*sinx);
	  this->f[j][1] += (this->f_hat[k_S][1]*cosx - this->f_hat[k_S][0]*sinx);
	}
    }
} /* void sparse_ndft_trafo */

void sparse_ndft_trafo_3d(sparse_nfft_plan_3d *this)
{
  int j,k_S,k0,k1,k2;
  double omega,cosx,sinx;
  int N=int_2_pow(this->J+2);
  int k_L;

  memset(this->f,0,this->M*sizeof(fftw_complex));

  for(k_S=0;k_S<this->N_S;k_S++)
    {
      k_L=this->index_sparse_to_full[k_S];

      k0=k_L/(N*N);
      k1=(k_L/N)%N;
      k2=k_L%N;
      
      for(j=0;j<this->M;j++)
	{
	  omega =
	    ((double)(k0 - N/2)) * this->act_nfft_plan->x[3 * j + 0] + 
	    ((double)(k1 - N/2)) * this->act_nfft_plan->x[3 * j + 1] +
	    ((double)(k2 - N/2)) * this->act_nfft_plan->x[3 * j + 2];
	  cosx = cos(2*PI*omega); sinx = sin(2*PI*omega);

	  this->f[j][0] += (this->f_hat[k_S][0]*cosx + this->f_hat[k_S][1]*sinx);
	  this->f[j][1] += (this->f_hat[k_S][1]*cosx - this->f_hat[k_S][0]*sinx);
	}
    }
} /* void sparse_ndft_trafo */


void sparse_nfft_trafo(sparse_nfft_plan *this)
{
  int r,rr,j;
  double temp,ctx,cty,stx,sty;

  int M=this->M;
  int J=this->J;

  /* center */
  this->center_nfft_plan->f_hat=this->f_hat+4*((J+1)/2+1)*int_2_pow(J);
  
  if (this->center_nfft_plan->N[0]<=this->center_nfft_plan->m) 
    ndft_trafo(this->center_nfft_plan);
  else
    nfft_trafo(this->center_nfft_plan);

  for (j=0; j<M; j++) 
    {
      this->f[j][0]= this->center_nfft_plan->f[j][0];
      this->f[j][1]= this->center_nfft_plan->f[j][1];
    }

  for(rr=0;rr<=(J+1)/2;rr++)
    {
      r=MIN(rr,J-rr);
      this->act_nfft_plan->my_fftw_plan1 = this->set_fftw_plan1[r];
      this->act_nfft_plan->N[0]=int_2_pow(r); this->act_nfft_plan->n[0]=this->sigma*this->act_nfft_plan->N[0];
      this->act_nfft_plan->N[1]=int_2_pow(J-r); this->act_nfft_plan->n[1]=this->sigma*this->act_nfft_plan->N[1];

      /*printf("%d x %d\n",this->act_nfft_plan->N[0],this->act_nfft_plan->N[1]);*/

      temp=-3.0*PI*int_2_pow(J-rr);

      /* right */      
      this->act_nfft_plan->f_hat=this->f_hat+(4*rr+0)*int_2_pow(J);
      
      if(r<rr)
	SWAP(this->act_nfft_plan->x,this->x_transposed);

      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  ndft_trafo(this->act_nfft_plan);
	else
	  short_nfft_trafo(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if(r<rr)
	SWAP(this->act_nfft_plan->x,this->x_transposed);
      
      for (j=0; j<M; j++)
	{
	  cty=cos(temp*this->act_nfft_plan->x[2*j+1]); sty=sin(temp*this->act_nfft_plan->x[2*j+1]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*cty - this->act_nfft_plan->f[j][1]*sty);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*cty + this->act_nfft_plan->f[j][0]*sty);
	}

      /* top */
      this->act_nfft_plan->f_hat=this->f_hat+(4*rr+1)*int_2_pow(J);
      
      if((r==rr)&&(J-rr!=rr))
	SWAP(this->act_nfft_plan->x,this->x_transposed);

      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  ndft_trafo(this->act_nfft_plan);
	else
	  short_nfft_trafo(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
      else
	nfft_trafo(this->act_nfft_plan);

      if((r==rr)&&(J-rr!=rr))
	SWAP(this->act_nfft_plan->x,this->x_transposed);

      for (j=0; j<M; j++)
	{
	  ctx=cos(temp*this->act_nfft_plan->x[2*j+0]); stx=sin(temp*this->act_nfft_plan->x[2*j+0]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ctx - this->act_nfft_plan->f[j][1]*stx);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ctx + this->act_nfft_plan->f[j][0]*stx);
	}

      /* left */
      this->act_nfft_plan->f_hat=this->f_hat+(4*rr+2)*int_2_pow(J);
      
      if(r<rr)
	SWAP(this->act_nfft_plan->x,this->x_transposed);

      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  ndft_trafo(this->act_nfft_plan);
	else
	  short_nfft_trafo(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
      else
	nfft_trafo(this->act_nfft_plan);

      if(r<rr)
	SWAP(this->act_nfft_plan->x,this->x_transposed);

      for (j=0; j<M; j++)
	{
	  cty=cos(temp*this->act_nfft_plan->x[2*j+1]); sty=sin(temp*this->act_nfft_plan->x[2*j+1]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*cty + this->act_nfft_plan->f[j][1]*sty);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*cty - this->act_nfft_plan->f[j][0]*sty);
	}

      /* bottom */
      this->act_nfft_plan->f_hat=this->f_hat+(4*rr+3)*int_2_pow(J);
      
      if((r==rr)&&(J-rr!=rr))
	SWAP(this->act_nfft_plan->x,this->x_transposed);

      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  ndft_trafo(this->act_nfft_plan);
	else
	  short_nfft_trafo(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if((r==rr)&&(J-rr!=rr))
	SWAP(this->act_nfft_plan->x,this->x_transposed);

      for (j=0; j<M; j++)
	{
	  ctx=cos(temp*this->act_nfft_plan->x[2*j+0]); stx=sin(temp*this->act_nfft_plan->x[2*j+0]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ctx + this->act_nfft_plan->f[j][1]*stx);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ctx - this->act_nfft_plan->f[j][0]*stx);
	}
    } /* for(rr) */
} /* void sparse_nfft_trafo */


void sparse_nfft_trafo_3d(sparse_nfft_plan_3d *this)
{
  int r,rr,j;
  double temp,ct0,ct1,ct2,st0,st1,st2;
  int sum_N_B_less_r,N_B_r,a,b;

  int M=this->M;
  int J=this->J;

  /* center */
  this->center_nfft_plan->f_hat=this->f_hat+6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1);
  
  if (this->center_nfft_plan->N[0]<=this->center_nfft_plan->m)
    ndft_trafo(this->center_nfft_plan);
  else
    nfft_trafo(this->center_nfft_plan);

  for (j=0; j<M; j++) 
    {
      this->f[j][0]= this->center_nfft_plan->f[j][0];
      this->f[j][1]= this->center_nfft_plan->f[j][1];
    }

  sum_N_B_less_r=0;
  for(rr=0;rr<=(J+1)/2;rr++)
    {
      a=int_2_pow(J-rr);
      b=int_2_pow(rr);

      N_B_r=a*b*b;

      r=MIN(rr,J-rr);
      this->act_nfft_plan->my_fftw_plan1 = this->set_fftw_plan1[rr];

      this->act_nfft_plan->N[0]=int_2_pow(r);
      if(a<b)
	this->act_nfft_plan->N[1]=int_2_pow(J-r);
      else
	this->act_nfft_plan->N[1]=int_2_pow(r);
      this->act_nfft_plan->N[2]=int_2_pow(J-r);

      /*printf("\n\n%d x %d x %d:\t",this->act_nfft_plan->N[0],this->act_nfft_plan->N[1],this->act_nfft_plan->N[2]); fflush(stdout);*/

      this->act_nfft_plan->N_L=this->act_nfft_plan->N[0]*this->act_nfft_plan->N[1]*this->act_nfft_plan->N[2];
      this->act_nfft_plan->n[0]=this->sigma*this->act_nfft_plan->N[0];
      this->act_nfft_plan->n[1]=this->sigma*this->act_nfft_plan->N[1];
      this->act_nfft_plan->n[2]=this->sigma*this->act_nfft_plan->N[2];
      this->act_nfft_plan->n_L=this->act_nfft_plan->n[0]*this->act_nfft_plan->n[1]*this->act_nfft_plan->n[2];

      /* only for right - rear - top */
      if((J==0)||((J==1)&&(rr==1)))
	temp=-2.0*PI;
      else
	temp=-3.0*PI*int_2_pow(J-rr);
	
      /* right */      
      this->act_nfft_plan->f_hat=this->f_hat + sum_N_B_less_r + N_B_r*0;

      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_120);

      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  if(this->act_nfft_plan->N[2]<=this->act_nfft_plan->m)
	    ndft_trafo(this->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(this->act_nfft_plan,&(this->set_nfft_plan_2d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_120);
      
      for (j=0; j<M; j++)
	{
	  ct0=cos(temp*this->act_nfft_plan->x[3*j+0]); st0=sin(temp*this->act_nfft_plan->x[3*j+0]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ct0 - this->act_nfft_plan->f[j][1]*st0);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ct0 + this->act_nfft_plan->f[j][0]*st0);
	}

      /* rear */      
      this->act_nfft_plan->f_hat=this->f_hat + sum_N_B_less_r + N_B_r*1;

      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_021);
      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_102);
      
      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  if(this->act_nfft_plan->N[2]<=this->act_nfft_plan->m)
	    ndft_trafo(this->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(this->act_nfft_plan,&(this->set_nfft_plan_2d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_021);
      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_102);
      
      for (j=0; j<M; j++)
	{
	  ct1=cos(temp*this->act_nfft_plan->x[3*j+1]); st1=sin(temp*this->act_nfft_plan->x[3*j+1]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ct1 - this->act_nfft_plan->f[j][1]*st1);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ct1 + this->act_nfft_plan->f[j][0]*st1);
	}

      /* top */      
      this->act_nfft_plan->f_hat=this->f_hat + sum_N_B_less_r + N_B_r*2;

      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_201);
      
      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  if(this->act_nfft_plan->N[2]<=this->act_nfft_plan->m)
	    ndft_trafo(this->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(this->act_nfft_plan,&(this->set_nfft_plan_2d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_201);
      
      for (j=0; j<M; j++)
	{
	  ct2=cos(temp*this->act_nfft_plan->x[3*j+2]); st2=sin(temp*this->act_nfft_plan->x[3*j+2]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ct2 - this->act_nfft_plan->f[j][1]*st2);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ct2 + this->act_nfft_plan->f[j][0]*st2);
	}

      /* only for left - front - bottom */
      if((J==0)||((J==1)&&(rr==1)))
	temp=-4.0*PI;
      else
	temp=-3.0*PI*int_2_pow(J-rr);

      /* left */
      this->act_nfft_plan->f_hat=this->f_hat + sum_N_B_less_r + N_B_r*3;

      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_120);

      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  if(this->act_nfft_plan->N[2]<=this->act_nfft_plan->m)
	    ndft_trafo(this->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(this->act_nfft_plan,&(this->set_nfft_plan_2d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_120);
      
      for (j=0; j<M; j++)
	{
	  ct0=cos(temp*this->act_nfft_plan->x[3*j+0]); st0=sin(temp*this->act_nfft_plan->x[3*j+0]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ct0 + this->act_nfft_plan->f[j][1]*st0);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ct0 - this->act_nfft_plan->f[j][0]*st0);
	}

      /* front */      
      this->act_nfft_plan->f_hat=this->f_hat + sum_N_B_less_r + N_B_r*4;

      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_021);
      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_102);
      
      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  if(this->act_nfft_plan->N[2]<=this->act_nfft_plan->m)
	    ndft_trafo(this->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(this->act_nfft_plan,&(this->set_nfft_plan_2d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if(a>b)
	SWAP(this->act_nfft_plan->x,this->x_021);
      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_102);
      
      for (j=0; j<M; j++)
	{
	  ct1=cos(temp*this->act_nfft_plan->x[3*j+1]); st1=sin(temp*this->act_nfft_plan->x[3*j+1]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ct1 + this->act_nfft_plan->f[j][1]*st1);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ct1 - this->act_nfft_plan->f[j][0]*st1);
	}

      /* bottom */      
      this->act_nfft_plan->f_hat=this->f_hat + sum_N_B_less_r + N_B_r*5;

      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_201);
      
      if(this->act_nfft_plan->N[0]<=this->act_nfft_plan->m)
	if(this->act_nfft_plan->N[1]<=this->act_nfft_plan->m)
	  if(this->act_nfft_plan->N[2]<=this->act_nfft_plan->m)
	    ndft_trafo(this->act_nfft_plan);
	  else
	    short_nfft_trafo_3d_1(this->act_nfft_plan,&(this->set_nfft_plan_1d[r]));
	else
	  short_nfft_trafo_3d_2(this->act_nfft_plan,&(this->set_nfft_plan_2d[r]));
      else
	nfft_trafo(this->act_nfft_plan);
      
      if(a<b)
	SWAP(this->act_nfft_plan->x,this->x_201);
      
      for (j=0; j<M; j++)
	{
	  ct2=cos(temp*this->act_nfft_plan->x[3*j+2]); st2=sin(temp*this->act_nfft_plan->x[3*j+2]);

	  this->f[j][0]+= (this->act_nfft_plan->f[j][0]*ct2 + this->act_nfft_plan->f[j][1]*st2);
	  this->f[j][1]+= (this->act_nfft_plan->f[j][1]*ct2 - this->act_nfft_plan->f[j][0]*st2);
	}

      sum_N_B_less_r+=6*N_B_r; 
    } /* for(rr) */
} /* void sparse_nfft_trafo_3d */

/*========================================================*/
/* J >1, no precomputation at all!! */
void sparse_nfft_init(sparse_nfft_plan *this, int J, int M, int m, unsigned snfft_flags)
{
  int r;
  int N[2];
  int n[2];

  this->snfft_flags=snfft_flags;
  this->sigma=2;
  this->J=J;
  this->M=M;
  this->N_S=(J+4)*int_2_pow(J+1);
  
  /* memory allocation */
  this->f = (fftw_complex *)fftw_malloc(M*sizeof(fftw_complex));
  this->f_hat = (fftw_complex *)fftw_malloc(this->N_S*sizeof(fftw_complex));
  this->x_transposed= (double*)fftw_malloc(2*M*sizeof(double));

  this->act_nfft_plan = (nfft_plan*)fftw_malloc(sizeof(nfft_plan));
  this->center_nfft_plan = (nfft_plan*)fftw_malloc(sizeof(nfft_plan));

  this->set_fftw_plan1=(fftw_plan*) fftw_malloc((J/2+1)*sizeof(fftw_plan));

  this->set_nfft_plan_1d = (nfft_plan*) fftw_malloc((ld(m)+1)*(sizeof(nfft_plan)));
  
  /* planning the small nffts */
  /* r=0 */
  N[0]=1;            n[0]=this->sigma*N[0];
  N[1]=int_2_pow(J); n[1]=this->sigma*N[1];
  
  nfft_init_specific(this->act_nfft_plan,2,N,M,n,m, PRE_LIN_PSI| MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
  
  if(this->act_nfft_plan->nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(this->act_nfft_plan,100000);

  this->set_fftw_plan1[0]=this->act_nfft_plan->my_fftw_plan1;

  for(r=1;r<=J/2;r++)
    {
      N[0]=int_2_pow(r);   n[0]=this->sigma*N[0];
      N[1]=int_2_pow(J-r); n[1]=this->sigma*N[1];
      this->set_fftw_plan1[r] = 
	fftw_plan_dft(2, n, this->act_nfft_plan->g1, this->act_nfft_plan->g2,
		      FFTW_FORWARD, this->act_nfft_plan->fftw_flags);
    }

  /* planning the 1d nffts */
  for(r=0;r<=ld(m);r++)
    {
      N[0]=int_2_pow(J-r); n[0]=this->sigma*N[0]; /* ==N[1] of the 2 dimensional plan */

      nfft_init_specific(&(this->set_nfft_plan_1d[r]),1,N,M,n,m, MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
      this->set_nfft_plan_1d[r].nfft_flags = this->set_nfft_plan_1d[r].nfft_flags | PRE_LIN_PSI;
      this->set_nfft_plan_1d[r].K=this->act_nfft_plan->K;
      this->set_nfft_plan_1d[r].psi=this->act_nfft_plan->psi;
    }

  /* center plan */
  /* J/2 == floor(((double)J) / 2.0) */
  N[0]=int_2_pow(J/2+1); n[0]=this->sigma*N[0];
  N[1]=int_2_pow(J/2+1); n[1]=this->sigma*N[1];
  nfft_init_specific(this->center_nfft_plan,2,N,M,n, m, MALLOC_F| FFTW_INIT, FFTW_MEASURE);
  this->center_nfft_plan->x= this->act_nfft_plan->x;
  this->center_nfft_plan->nfft_flags = this->center_nfft_plan->nfft_flags | PRE_LIN_PSI;
  this->center_nfft_plan->K=this->act_nfft_plan->K;
  this->center_nfft_plan->psi=this->act_nfft_plan->psi;

  if(this->snfft_flags & SNDFT)
    {
      this->index_sparse_to_full=(int*)fftw_malloc(this->N_S*sizeof(int));
      init_index_sparse_to_full(this);
    }
}

/*========================================================*/
/* J >1, no precomputation at all!! */
void sparse_nfft_init_3d(sparse_nfft_plan_3d *this, int J, int M, int m, unsigned snfft_flags)
{
  int r,rr,a,b;
  int N[3];
  int n[3];

  this->snfft_flags=snfft_flags;
  this->sigma=2;
  this->J=J;
  this->M=M;
  this->N_S=6*int_2_pow(J)*(int_2_pow((J+1)/2+1)-1)+int_2_pow(3*(J/2+1));
  
  /* memory allocation */
  this->f =     (fftw_complex *)fftw_malloc(M*sizeof(fftw_complex));
  this->f_hat = (fftw_complex *)fftw_malloc(this->N_S*sizeof(fftw_complex));

  this->x_102= (double*)fftw_malloc(3*M*sizeof(double));
  this->x_201= (double*)fftw_malloc(3*M*sizeof(double));
  this->x_120= (double*)fftw_malloc(3*M*sizeof(double));
  this->x_021= (double*)fftw_malloc(3*M*sizeof(double));

  this->act_nfft_plan = (nfft_plan*)fftw_malloc(sizeof(nfft_plan));
  this->center_nfft_plan = (nfft_plan*)fftw_malloc(sizeof(nfft_plan));

  this->set_fftw_plan1=(fftw_plan*) fftw_malloc(((J+1)/2+1)*sizeof(fftw_plan));

  this->set_nfft_plan_1d = (nfft_plan*) fftw_malloc((ld(m)+1)*(sizeof(nfft_plan)));
  this->set_nfft_plan_2d = (nfft_plan*) fftw_malloc((ld(m)+1)*(sizeof(nfft_plan)));
  
  /* planning the small nffts */
  /* r=0 */
  N[0]=1;            n[0]=this->sigma*N[0];
  N[1]=1;            n[1]=this->sigma*N[1];
  N[2]=int_2_pow(J); n[2]=this->sigma*N[2];
  
  nfft_init_specific(this->act_nfft_plan,3,N,M,n,m, PRE_LIN_PSI| MALLOC_X| MALLOC_F, FFTW_MEASURE);

  if(this->act_nfft_plan->nfft_flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(this->act_nfft_plan,100000);

  /* malloc g1, g2 for maximal size */
  this->act_nfft_plan->g1 = fftw_malloc(this->sigma*this->sigma*this->sigma*int_2_pow(J+(J+1)/2)*sizeof(fftw_complex));
  this->act_nfft_plan->g2 = fftw_malloc(this->sigma*this->sigma*this->sigma*int_2_pow(J+(J+1)/2)*sizeof(fftw_complex));

  this->act_nfft_plan->my_fftw_plan1 =
    fftw_plan_dft(3, n, this->act_nfft_plan->g1, this->act_nfft_plan->g2,
		  FFTW_FORWARD, this->act_nfft_plan->fftw_flags);

  this->set_fftw_plan1[0]=this->act_nfft_plan->my_fftw_plan1;

  for(rr=1;rr<=(J+1)/2;rr++)
    {
      a=int_2_pow(J-rr);
      b=int_2_pow(rr);

      r=MIN(rr,J-rr);

      n[0]=this->sigma*int_2_pow(r);
      if(a<b)
	n[1]=this->sigma*int_2_pow(J-r);
      else
	n[1]=this->sigma*int_2_pow(r);
      n[2]=this->sigma*int_2_pow(J-r);

      this->set_fftw_plan1[rr] =
	fftw_plan_dft(3, n, this->act_nfft_plan->g1, this->act_nfft_plan->g2,
		      FFTW_FORWARD, this->act_nfft_plan->fftw_flags);
    }

  /* planning the 1d nffts */
  for(r=0;r<=ld(m);r++)
    {
      N[0]=int_2_pow(J-r); n[0]=this->sigma*N[0];
      N[1]=int_2_pow(J-r); n[1]=this->sigma*N[1];

      if(N[0]>m)
	{
	  nfft_init_specific(&(this->set_nfft_plan_1d[r]),1,N,M,n,m, MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
	  this->set_nfft_plan_1d[r].nfft_flags = this->set_nfft_plan_1d[r].nfft_flags | PRE_LIN_PSI;
	  this->set_nfft_plan_1d[r].K=this->act_nfft_plan->K;
	  this->set_nfft_plan_1d[r].psi=this->act_nfft_plan->psi;
	  nfft_init_specific(&(this->set_nfft_plan_2d[r]),2,N,M,n,m, MALLOC_X| MALLOC_F| FFTW_INIT, FFTW_MEASURE);
	  this->set_nfft_plan_2d[r].nfft_flags = this->set_nfft_plan_2d[r].nfft_flags | PRE_LIN_PSI;
	  this->set_nfft_plan_2d[r].K=this->act_nfft_plan->K;
	  this->set_nfft_plan_2d[r].psi=this->act_nfft_plan->psi;
	}
    }

  /* center plan */
  /* J/2 == floor(((double)J) / 2.0) */
  N[0]=int_2_pow(J/2+1); n[0]=this->sigma*N[0];
  N[1]=int_2_pow(J/2+1); n[1]=this->sigma*N[1];
  N[2]=int_2_pow(J/2+1); n[2]=this->sigma*N[2];
  nfft_init_specific(this->center_nfft_plan,3,N,M,n, m, MALLOC_F| FFTW_INIT, FFTW_MEASURE);
  this->center_nfft_plan->x= this->act_nfft_plan->x;
  this->center_nfft_plan->nfft_flags = this->center_nfft_plan->nfft_flags | PRE_LIN_PSI;
  this->center_nfft_plan->K=this->act_nfft_plan->K;
  this->center_nfft_plan->psi=this->act_nfft_plan->psi;

  if(this->snfft_flags & SNDFT)
    {
      this->index_sparse_to_full=(int*)fftw_malloc(this->N_S*sizeof(int));
      init_index_sparse_to_full_3d(this);
    }
}

void sparse_nfft_finalize(sparse_nfft_plan *this)
{
  int r;

  if(this->snfft_flags & SNDFT)
    fftw_free(this->index_sparse_to_full);

  /* center plan */
  this->center_nfft_plan->nfft_flags = this->center_nfft_plan->nfft_flags ^ PRE_LIN_PSI;
  nfft_finalize(this->center_nfft_plan);

  /* the 1d nffts */
  for(r=0;r<=ld(this->act_nfft_plan->m);r++)
    {
      this->set_nfft_plan_1d[r].nfft_flags = this->set_nfft_plan_1d[r].nfft_flags ^ PRE_LIN_PSI;
      nfft_finalize(&(this->set_nfft_plan_1d[r]));
    }
  
  /* finalize the small nffts */
  this->act_nfft_plan->my_fftw_plan1=this->set_fftw_plan1[0];

  for(r=1;r<=this->J/2;r++)
    fftw_destroy_plan(this->set_fftw_plan1[r]);

  /* r=0 */
  nfft_finalize(this->act_nfft_plan);

  fftw_free(this->set_nfft_plan_1d);

  fftw_free(this->set_fftw_plan1);

  fftw_free(this->x_transposed);

  fftw_free(this->f_hat);
  fftw_free(this->f);
}

void sparse_nfft_finalize_3d(sparse_nfft_plan_3d *this)
{
  int r;

  if(this->snfft_flags & SNDFT)
    fftw_free(this->index_sparse_to_full);

  /* center plan */
  this->center_nfft_plan->nfft_flags = this->center_nfft_plan->nfft_flags ^ PRE_LIN_PSI;
  nfft_finalize(this->center_nfft_plan);

  /* the 1d and 2d nffts */
  for(r=0;r<=ld(this->act_nfft_plan->m);r++)
    {
      if(int_2_pow(this->J-r)>this->act_nfft_plan->m)
	{
	  this->set_nfft_plan_2d[r].nfft_flags = this->set_nfft_plan_2d[r].nfft_flags ^ PRE_LIN_PSI;
	  nfft_finalize(&(this->set_nfft_plan_2d[r]));
	  this->set_nfft_plan_1d[r].nfft_flags = this->set_nfft_plan_1d[r].nfft_flags ^ PRE_LIN_PSI;
	  nfft_finalize(&(this->set_nfft_plan_1d[r]));
	}
    }
  
  /* finalize the small nffts */
  this->act_nfft_plan->my_fftw_plan1=this->set_fftw_plan1[0];

  for(r=1;r<=(this->J+1)/2;r++)
    fftw_destroy_plan(this->set_fftw_plan1[r]);

  /* r=0 */
  nfft_finalize(this->act_nfft_plan);

  fftw_free(this->set_nfft_plan_1d);
  fftw_free(this->set_nfft_plan_2d);

  fftw_free(this->set_fftw_plan1);

  fftw_free(this->x_102);
  fftw_free(this->x_201);
  fftw_free(this->x_120);
  fftw_free(this->x_021);

  fftw_free(this->f_hat);
  fftw_free(this->f);
}
