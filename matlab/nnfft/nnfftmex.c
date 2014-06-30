/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: nnfftmex.c 4086 2014-06-23 19:09:50Z sukunis $ */
#include "config.h"

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "nfft3.h"
#include "infft.h"
#include "imex.h"

#ifdef HAVE_MEXVERSION_C
  #include "mexversion.c"
#endif

#define PLANS_MAX 100 /* maximum number of plans */
#define CMD_LEN_MAX 20 /* maximum length of command argument */

/* global flags */
#define NNFFT_MEX_FIRST_CALL (1U << 0)
unsigned short gflags = NNFFT_MEX_FIRST_CALL;

nnfft_plan* plans[PLANS_MAX]; /* plans */
char cmd[CMD_LEN_MAX];

static inline void get_nm(const mxArray *prhs[], int *n, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"nnfft: Input argument N must be a scalar.");
  DM(if ((t < 0) || (t%2!=0))
    mexErrMsgTxt("nnfft: Input argument N must be non-negative and multiple of two.");)
  *n = t;
  t = nfft_mex_get_int(prhs[2],"nnfft: Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("nnfft: Input argument M must be positive.");)
  *m = t;
}




static inline void check_nargs(const int nrhs, const int n, const char* errmsg)
{
  DM(if (nrhs != n)
    mexErrMsgTxt(errmsg);)
}

static inline void check_plan(int i)
{
  DM(if (i < 0 || i >= PLANS_MAX)
    mexErrMsgTxt("Invalid plan");)
  DM(if (plans[i] == 0)
    mexErrMsgTxt("Plan was not initialized or has already been finalized");)
}

static inline int mkplan(void)
{
  mexLock();
  int i = 0;
  while (i < PLANS_MAX && plans[i] != 0) i++;
  if (i == PLANS_MAX)
    mexErrMsgTxt("nnfft: Too many plans already allocated.");
  plans[i] = nfft_malloc(sizeof(nnfft_plan));
  return i;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  mexUnlock();

  if (!(gflags & NNFFT_MEX_FIRST_CALL))
  {
    for (i = 0; i < PLANS_MAX; i++)
      if (plans[i])
      {
        nnfft_finalize(plans[i]);
        plans[i] = 0;
      }
    gflags |= NNFFT_MEX_FIRST_CALL;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (gflags & NNFFT_MEX_FIRST_CALL)
  {
    /* Force Matlab to load libfftw3. There is at least one version of Matlab
     * which otherwise crashes upon invocation of this mex function. */
    mexEvalString("fft([1,2,3,4]);");

    nfft_mex_install_mem_hooks();

    /* plan pointers to zeros */
    {
      int i;
      for (i = 0; i < PLANS_MAX; i++)
        plans[i] = 0;
    }

    mexAtExit(cleanup);
    gflags &= ~NNFFT_MEX_FIRST_CALL;
  }

  /* command string */
  DM(if (nrhs == 0)
    mexErrMsgTxt("At least one input required.");)

  DM(if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument must be a string.");)

  if (mxGetString(prhs[0], cmd, CMD_LEN_MAX))
    mexErrMsgTxt("Could not get command string.");

  if (strcmp(cmd,"init_1d") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for init.");
    {
      int i;
      int n, m;
      get_nm(prhs,&n,&m);
      disp("ok");
      i = mkplan();
      nnfft_init_1d(plans[i],n,m);
      plhs[0] = mxCreateDoubleScalar((double)i);
    }
    return;
  }
  else if (strcmp(cmd,"init_2d") == 0)
  {
	  mexErrMsgTxt("No implementation available.");
//    check_nargs(nrhs,4,"Wrong number of arguments for init.");
//    {
//      int i;
//      int n1, n2, m;
//      get_n1n2m(prhs,&n1,&n2,&m);
//      i = mkplan();
//      nfft_init_2d(plans[i],n1,n2,m);
//      plhs[0] = mxCreateDoubleScalar((double)i);
//    }
//    return;
  }
  else if (strcmp(cmd,"init_3d") == 0)
  {
	  mexErrMsgTxt("No implementation available.");
//    check_nargs(nrhs,5,"Wrong number of arguments for init.");
//    {
//      int i;
//      int n1, n2, n3, m;
//      get_n1n2n3m(prhs,&n1,&n2,&n3,&m);
//      i = mkplan();
//      nfft_init_3d(plans[i],n1,n2,n3,m);
//      plhs[0] = mxCreateDoubleScalar((double)i);
//    }
//    return;
  }
  else if (strcmp(cmd,"init_guru") == 0)
  {
	  mexErrMsgTxt("No implementation available.");
//    /** NO ERROR HANDLING !!*/
//    int i;
//    const int d = mxGetScalar(mxGetCell(prhs[1], 0));
//
//    int N[d],n[d],m,M;
//    unsigned int f1,f2;
//
//    DM(if ((d < 1) || (d>4))
//	 mexErrMsgTxt("nfft: Input argument d must be positive and smaller than 5.");)
//
//    get_guru(prhs,d,N,&M,n,&m,&f1,&f2);
//    i = mkplan();
//    nfft_init_guru(plans[i],d,N,M,n,m,
//		   f1 | MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT,
//		   f2);
//
//    plhs[0] = mxCreateDoubleScalar((double)i);
//
//    return;
  }
  else if (strcmp(cmd,"precompute_psi") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for precompute_one_psi.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      nnfft_precompute_psi(plans[i]);// hier die gewuenschte fkt fuer vorberechnung eintragen
    }
    return;
  }
  else if (strcmp(cmd,"trafo") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for trafo.");
    {
      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
      check_plan(i);
      nnfft_trafo(plans[i]);
    }
    return;
  }
  else if (strcmp(cmd,"adjoint") == 0)
  {
	  mexErrMsgTxt("No implementation available.");
//    check_nargs(nrhs,2,"Wrong number of arguments for adjoint.");
//    {
//      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
//      check_plan(i);
//      nfft_adjoint(plans[i]);
//    }
//    return;
  }
  else if (strcmp(cmd,"finalize") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for finalize.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      check_plan(i);
      nnfft_finalize(plans[i]);
      nfft_free(plans[i]);
      plans[i] = 0;
    }
    return;
  }
  else if (strcmp(cmd,"trafo_direct") == 0)
  {
	  mexErrMsgTxt("No implementation available.");
//    check_nargs(nrhs,2,"Wrong number of arguments for trafo direct.");
//    {
//      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
//      check_plan(i);
//      nnfft_trafo_direct(plans[i]);
//    }
//    return;
  }
  else if (strcmp(cmd,"adjoint_direct") == 0)
  {
	  mexErrMsgTxt("No implementation available.");
//    check_nargs(nrhs,2,"Wrong number of arguments for adjoint direct.");
//    {
//      int i = nfft_mex_get_int(prhs[1],"nfft: Input argument plan must be a scalar.");
//      check_plan(i);
//      nfft_adjoint_direct(plans[i]);
//    }
//    return;
  }
  else if (strcmp(cmd,"get_x") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_x.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      int m, d;
      check_plan(i);
      m = plans[i]->M_total;
      d = plans[i]->d;
      plhs[0] = mxCreateDoubleMatrix((unsigned int)d, (unsigned int)m, mxREAL);
      {
        double *x = mxGetPr(plhs[0]);
        int j,t;
        for (j = 0; j < m; j++)
	  for (t = 0; t < d; t++)
	    x[d*j+t] = plans[i]->x[d*j+t];
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_f.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      int m;
      check_plan(i);
      m = plans[i]->M_total;
      plhs[0] = mxCreateDoubleMatrix((unsigned int)m, 1, mxCOMPLEX);
      {
        double *fr = mxGetPr(plhs[0]), *fi = mxGetPi(plhs[0]);
        int j;
        for (j = 0; j < m; j++)
        {
          fr[j] = creal(plans[i]->f[j]);
          fi[j] = cimag(plans[i]->f[j]);
        }
      }
    }
    return;
  }
  else if (strcmp(cmd,"get_f_hat") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for get_f_hat.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      int n;
      check_plan(i);
      n = plans[i]->N_total;
      plhs[0] = mxCreateDoubleMatrix((unsigned int)n, 1, mxCOMPLEX);
      {
        double *f_hatr = mxGetPr(plhs[0]), *f_hati = mxGetPi(plhs[0]);
        int k;
        for (k = 0; k < n; k++)
        {
          f_hatr[k] = creal(plans[i]->f_hat[k]);
          f_hati[k] = cimag(plans[i]->f_hat[k]);
        }
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_x") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_x.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      int m, d;
      check_plan(i);
      m = plans[i]->M_total;
      d = plans[i]->d;
      DM(if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) > 2)
        mexErrMsgTxt("Input argument x must be a d x M double array");)
      DM(if (mxGetM(prhs[2]) != (unsigned int)d || mxGetN(prhs[2]) != (unsigned int)m)
        mexErrMsgTxt("Input argument x must have correct size (d x M).");)
      {
        double *x = mxGetPr(prhs[2]);
        int j,t;
        for (j = 0; j < m; j++)
	  for (t = 0; t < d; t++)
	    plans[i]->x[d*j+t] = x[d*j+t];
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_f") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_f.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      int m;
      check_plan(i);
      m = plans[i]->M_total;
      DM(if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != (unsigned int)m)
        mexErrMsgTxt("Input argument f must have correct size.");)
      {
        double *fr = mxGetPr(prhs[2]), *fi = mxGetPi(prhs[2]);
        int j;
        if (fi)
          for (j = 0; j < m; j++)
            plans[i]->f[j] = fr[j] + _Complex_I*fi[j];
        else
          for (j = 0; j < m; j++)
            plans[i]->f[j] = fr[j];
      }
    }
    return;
  }
  else if (strcmp(cmd,"set_f_hat") == 0)
  {
    check_nargs(nrhs,3,"Wrong number of arguments for set_f_hat.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      int n;
      check_plan(i);
      n = plans[i]->N_total;
      DM(if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("Input argument f must be a double array");)
      DM(if (   mxGetM(prhs[2]) != (unsigned int)n || mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Input argument f must have correct size.");)
      {
        double *f_hatr = mxGetPr(prhs[2]), *f_hati = mxGetPi(prhs[2]);
        int k;
        if (f_hati)
          for (k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k] + _Complex_I*f_hati[k];
        else
          for (k = 0; k < n; k++)
            plans[i]->f_hat[k] = f_hatr[k];
      }
    }
    return;
  }
  else if (strcmp(cmd,"display") == 0)
  {
    check_nargs(nrhs,2,"Wrong number of arguments for set_f_hat_linear.");
    {
      int i = nfft_mex_get_int(prhs[1],"nnfft: Input argument plan must be a scalar.");
      mexPrintf("Plan %d\n",i);
      check_plan(i);
      mexPrintf("  pointer: %p\n",plans[i]);
      mexPrintf("        d: %d\n",plans[i]->d);
      mexPrintf("  N_total: %d\n",plans[i]->N_total);
      mexPrintf("  M_total: %d\n",plans[i]->M_total);
      mexPrintf("        x: %p\n",plans[i]->x);
      mexPrintf("        f: %p\n",plans[i]->f);
      mexPrintf("    f_hat: %p\n",plans[i]->f_hat);
      //mexPrintf("    flags: %d\n",plans[i]->flags);
    }
    return;
  }
  else if(strcmp(cmd,"get_num_threads") == 0)
  {
    int32_t nthreads = X(get_num_threads)();
    plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *((int32_t *)mxGetData(plhs[0])) = nthreads;

    return;
  }
  else
    mexErrMsgTxt("nnfft: Unknown command.\n");
}
