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

#include <stdlib.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"

/**
 * \defgroup applications_mri3d_construct_data_3d construct_data_3d
 * \ingroup applications_mri3d
 * \{
 */

static void construct(char * file, int N, int M, int Z)
{
  int j,k,l;                /* some variables */
  double real;
  nfft_plan my_plan;        /* plan for the three dimensional nfft  */
  FILE* fp,*fk;
  int my_N[3],my_n[3];      /* to init the nfft */


  /* initialise my_plan */
  //nfft_init_3d(&my_plan,Z,N,N,M);
  my_N[0]=Z; my_n[0]=ceil(Z*1.2);
  my_N[1]=N; my_n[1]=ceil(N*1.2);
  my_N[2]=N; my_n[2]=ceil(N*1.2);
  nfft_init_guru(&my_plan, 3, my_N, M, my_n, 6,
                      PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT,
                      FFTW_MEASURE);

  fp=fopen("knots.dat","r");

  for(j=0;j<M;j++)
    fscanf(fp,"%le %le %le",&my_plan.x[3*(j)+1],
      &my_plan.x[3*(j)+2],&my_plan.x[3*(j)+0]);

  fclose(fp);

  fp=fopen("input_f.dat","r");
  fk=fopen(file,"w");

  for(l=0;l<Z;l++) {
    for(j=0;j<N;j++)
    {
      for(k=0;k<N;k++)
      {
        //fscanf(fp,"%le ",&my_plan.f_hat[(N*N*(Z-l)+N*j+k+N*N*Z/2)%(N*N*Z)][0]);
        fscanf(fp,"%le ",&real);
        my_plan.f_hat[(N*N*l+N*j+k)] = real;
      }
    }
  }

    if(my_plan.flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);

    nfft_trafo(&my_plan);


    for(j=0;j<my_plan.M_total;j++)
      fprintf(fk,"%le %le %le %le %le\n",my_plan.x[3*j+1],
      my_plan.x[3*j+2],my_plan.x[3*j+0],creal(my_plan.f[j]),cimag(my_plan.f[j]));



  fclose(fk);
  fclose(fp);

  nfft_finalize(&my_plan);
}

int main(int argc, char **argv)
{
  if (argc <= 4) {
    printf("usage: ./construct_data FILENAME N M Z\n");
    return 1;
  }

  construct(argv[1], atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));

	return 1;
}
/* \} */
