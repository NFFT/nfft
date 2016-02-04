/*
 * Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
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
#define NFSOFT_MEX_FIRST_CALL (1U << 0)
unsigned short gflags = NFSOFT_MEX_FIRST_CALL;

nfsoft_plan* plans[PLANS_MAX]; /* plans */
int n_max = -1; /* maximum degree precomputed */
char cmd[CMD_LEN_MAX];

static inline void get_nm(const mxArray *prhs[], int *n, int *m)
{
  int t = nfft_mex_get_int(prhs[1],"nfsoft: Input argument N must be a scalar.");
  DM(if (t < 0)
    mexErrMsgTxt("nfsoft: Input argument N must be non-negative.");)
  *n = t;
  t = nfft_mex_get_int(prhs[2],"nfsoft: Input argument M must be a scalar.");
  DM(if (t < 1)
    mexErrMsgTxt("nfsoft: Input argument M must be positive.");)
  *m = t;
}

static inline void get_nmf(const mxArray *prhs[], int *n, int *m,
  unsigned int *f)
{
  get_nm(prhs,n,m);
  *f = nfft_mex_get_int(prhs[3],"nfsoft: Input argument flags must be a scalar.");
}

static inline void get_nmffc(const mxArray *prhs[], int *n, int *m,
  unsigned int *f, unsigned int *f2, int *c)
{
  get_nmf(prhs,n,m,f);
  *f2 = nfft_mex_get_int(prhs[4],"nfsoft: Input argument flags2 must be a scalar.");
  {
    int t = nfft_mex_get_int(prhs[5],"nfsoft: Input argument c must be a scalar.");
    DM(if (t < 1)
      mexErrMsgTxt("nfsoft: Input argument c must be positive.");)
    *c = t;
  }
}

static inline void check_nargs(const int nrhs, const int n, const char* errmsg)
{
  DM(if (nrhs != n)
    mexErrMsgTxt(errmsg);)
}

static inline int mkplan()
{
  int i = 0;
  while (i < PLANS_MAX && plans[i] != 0) i++;
  if (i == PLANS_MAX)
    mexErrMsgTxt("nfsoft: Too many plans already allocated.");
  plans[i] = nfft_malloc(sizeof(nfsoft_plan));
  return i;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  if (!(gflags & NFSOFT_MEX_FIRST_CALL))
  {
    for (i = 0; i < PLANS_MAX; i++)
      if (plans[i])
      {
        nfsoft_finalize(plans[i]);
        nfft_free(plans[i]);
        plans[i] = 0;
      }
    gflags |= NFSOFT_MEX_FIRST_CALL;
    n_max = -1;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (gflags & NFSOFT_MEX_FIRST_CALL)
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
    gflags &= ~NFSOFT_MEX_FIRST_CALL;
  }

  /* command string */
  DM(if (nrhs == 0)
    mexErrMsgTxt("At least one input required.");)

  DM(if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument must be a string.");)

  if (mxGetString(prhs[0], cmd, CMD_LEN_MAX))
    mexErrMsgTxt("Could not get command string.");

  if(strcmp(cmd,"get_num_threads") == 0)
  {
    int32_t nthreads = X(get_num_threads)();
    plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *((int32_t *)mxGetData(plhs[0])) = nthreads;

    return;
  }
  else
    mexErrMsgTxt("nfsoft: Unknown command.\n");
}

