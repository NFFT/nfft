#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)

int int_2_pow(int a)
{
  return (1U<< a);
}

int index_sparse_to_full(int J, int k)
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

	if (o==0)
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
	}
	else if (o==1)
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
	}
	else if (o==2)
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
	}
	else if (o==3)
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
	}
    }

    //printf("m1=%d, m2=%d\n",m1,m2);

    return(k1*N+k2);
}

int index_full_to_sparse(int J, int k)
{
    int N=int_2_pow(J+2);  /* number of full coeffs       */
    int N_B=int_2_pow(J);  /* number in each sparse block */

    int k1=k/N-N/2;        /* coordinates in the full grid */
    int k2=k%N-N/2;        /* k1: row, k2: column          */

    int r,a,b,o=-1;

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

int main(int argc, char *argv[])
{
    int J, k, j, N;

    /* evaluate parameters */
    if (argc!=3)  {
	printf("index_func J k \n");
	exit(-1);
    }
    
    J = atoi(argv[1]);
    k = atoi(argv[2]);
    N = int_2_pow(J+2);
    
    j=index_sparse_to_full(J,k);

    printf("sparse index %d bei J=%d ist full index %d, entspricht (%d,%d)\n",k,J,j,j/N-N/2,j%N-N/2);

    printf("sparse index = %d\n",index_full_to_sparse(J,index_sparse_to_full(J,k)));

    for (j=0; j<N; j++)
    {
	for (k=0; k<N; k++)
	{
            int temp=index_full_to_sparse(J,j*N+k);
	    if (temp==-1)
		printf(" ");
	    else
		printf("%d",temp%10);
	}
	printf("\n");
    }

    return 0;
}
