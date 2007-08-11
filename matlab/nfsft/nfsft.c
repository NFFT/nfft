/* $Id$
 *
 * Copyright (c) 2007 Jens Keiner
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

#include <math.h>
#include <string.h>
#include "nfft3.h"
#include "util.h"
#include "../imex.h"
/*----------------------------------------------------------------------------*/

#ifndef HAVE_MEXVERSION_C
  #include "mexversion.c"
#endif
/*----------------------------------------------------------------------------*/

#define PLANS_MAX 10             /**< Maximum number of plans.                */
#define CMD_LEN_MAX 20           /**< Maximum length of command argument.     */
/*----------------------------------------------------------------------------*/

unsigned short first_call = 1U;  /**< Flag for first call to mexFunction      */
unsigned short precomputed = 0U; /**< Flag for nfsft_precomputed invocation   */
nfsft_plan* plans[PLANS_MAX];    /**< Plan pointers array                     */
int n_max = -1;                  /**< Maximum degree after precomputation     */
/*----------------------------------------------------------------------------*/

/** Gets and checks degree and number of nodes. */
#define ARG_GET_NM(x,y) \
ARG_GET_NONNEG_INT(1,Second,x) \
ARG_GET_POS_INT(2,Third,y) \
if (!precomputed || (int)(x[0]) > n_max) \
  mexErrMsgTxt("Degree too high (execute precompute with high enough degree).");

/** Gets nfsft flags. */
#define ARG_GET_FLAGS(x) \
ARG_GET_NONNEG_INT(3,Fourth,x)

/** Gets nfft flags and cutoff. */
#define ARG_GET_FLAGS2(x,y) \
ARG_GET_NONNEG_INT(4,Fifth,x) \
ARG_GET_NONNEG_INT(5,Sixth,y)

/** Gets plan reference number. */
#define ARG_GET_PLAN(x) \
ARG_GET_NONNEG_INT(1,Second,x) \
if (x[0] >= (double)PLANS_MAX || plans[(int)(x[0])] == 0) \
  mexErrMsgTxt("Invalid plan number.");

/** Creates a new plan. */
#define MAKE_NEW_PLAN(x) \
x = 0; \
while (x < PLANS_MAX && plans[x] != 0) \
{ \
  x++; \
} \
if (x == PLANS_MAX) \
  mexErrMsgTxt("Too many plans already allocated"); \
plans[x] = nfft_malloc(sizeof(nfsft_plan));

/** Returns the plan reference number. */
#define RETURN_PLAN(x,y) \
plhs[0] = mxCreateDoubleScalar(mxREAL); \
x = mxGetPr(plhs[0]); \
x[0] = y;
/*----------------------------------------------------------------------------*/

void cleanup(void)
{
  int i;

  if (!first_call)
  {
    for (i = 0; i < PLANS_MAX; i++)
    {
      if (plans[i] != 0)
      {
        nfsft_finalize(plans[i]);
        nfft_free(plans[i]);
        plans[i] = 0;
      }
    }
    nfsft_forget();
    first_call = 1U;
    precomputed = 0U;
    n_max = -1;
  }
}
/*----------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int cmdl;                            /**< Length of string buffer           */
  char *cmd;                           /**< Command string                    */
  int i, j;                            /**< Loop variable                     */
  double *dp1, *dp2, *dp3, *dp4, *dp5; /**< Argument pointers                 */

  if (first_call)
  {
    install_mem_hooks();

    /* Initialize plan pointers with zeros. */
    for (i = 0; i < PLANS_MAX; i++)
    {
      plans[i] = 0;
    }

    mexAtExit(cleanup);
    first_call = 0U;
  }

  /* Check for command name. */
  if (nrhs == 0)
  {
    mexErrMsgTxt("At least one input required.");
  }

  /* Check if first argument is a string. */
  if (mxIsChar(prhs[0]) != 1)
  {
    mexErrMsgTxt("First argument must be a string.");
  }

  /* Command string length */
  cmdl = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

  /* Check for too long command string. */
  if (cmdl > CMD_LEN_MAX+1)
  {
    mexErrMsgTxt("Command argument too long.");
  }

  /* Allocate buffer for the command string. */
  cmd = mxCalloc(cmdl, sizeof(char));

  /* Copy the string data to buffe. */
  if (mxGetString(prhs[0], cmd, cmdl) != 0)
  {
    mexErrMsgTxt("Could not get command string.");
  }

  /* Command execution */
  if (strcmp(cmd,"init") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 3)
      mexErrMsgTxt("Two additional arguments required for init.");

    ARG_GET_NM(dp1,dp2)
    MAKE_NEW_PLAN(i)
    nfsft_init(plans[i],(int)(dp1[0]),(int)(dp2[0]));
    RETURN_PLAN(dp1,i)
  }
  else if (strcmp(cmd,"init_advanced") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 4)
      mexErrMsgTxt("Three additional arguments required for init_advanced.");

    ARG_GET_NM(dp1,dp2)
    ARG_GET_FLAGS(dp3)
    MAKE_NEW_PLAN(i)
    nfsft_init_advanced(plans[i],(int)(dp1[0]),(int)(dp2[0]),
      (unsigned int)(dp3[0]) | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
      NFSFT_MALLOC_F_HAT);
    RETURN_PLAN(dp1,i)
  }
  else if (strcmp(cmd,"init_guru") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 6)
      mexErrMsgTxt("Five additional arguments required for init_guru.");

    ARG_GET_NM(dp1,dp2)
    ARG_GET_FLAGS(dp3)
    ARG_GET_FLAGS2(dp4,dp5)
    MAKE_NEW_PLAN(i)
    nfsft_init_guru(plans[i],(int)(dp1[0]),(int)(dp2[0]),
      (unsigned int)(dp3[0]) | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
      NFSFT_MALLOC_F_HAT, (int)(dp4[0]),(int)(dp5[0]));
    RETURN_PLAN(dp1,i)
  }
  else if (strcmp(cmd,"precompute") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 5)
      mexErrMsgTxt("Four additional arguments required for precompute.");

    ARG_GET_NONNEG_INT(1,Second,dp1)
    ARG_GET_NONNEG_DOUBLE(2,Third,dp2)
    ARG_GET_NONNEG_INT(3,Fourth,dp3)
    ARG_GET_NONNEG_INT(4,Fifth,dp4)

    if (n_max < (int)(dp1[0]))
    {
      if (precomputed)
        nfsft_forget();

      nfsft_precompute((int)(dp1[0]),dp2[0],(unsigned int)(dp3[0]),
        (unsigned int)(dp4[0]));
      precomputed = 1U;
      n_max = (int)(dp1[0]);
    }
  }
  else if (strcmp(cmd,"forget") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 1)
      mexErrMsgTxt("No further arguments allowed for forget.");

    nfsft_forget();
  }
  else if (strcmp(cmd,"trafo") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for trafo.");

    ARG_GET_PLAN(dp1)
    nfsft_trafo(plans[(int)(dp1[0])]);
  }
  else if (strcmp(cmd,"adjoint") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for adjoint.");

    ARG_GET_PLAN(dp1)
    nfsft_adjoint(plans[(int)(dp1[0])]);
  }
  else if (strcmp(cmd,"finalize") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for finalize.");

    ARG_GET_PLAN(dp1)
    nfsft_finalize(plans[(int)(dp1[0])]);
    nfft_free(plans[(int)(dp1[0])]);
    plans[(int)(dp1[0])] = 0;
  }
  else if (strcmp(cmd,"precompute_x") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for precompute_x.");

    ARG_GET_PLAN(dp1)
    nfsft_precompute_x(plans[(int)(dp1[0])]);
  }
  else if (strcmp(cmd,"trafo_direct") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for trafo_direct.");

    ARG_GET_PLAN(dp1)
    ndsft_trafo(plans[(int)(dp1[0])]);
  }
  else if (strcmp(cmd,"adjoint_direct") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for adjoint_direct.");

    ARG_GET_PLAN(dp1)
    ndsft_adjoint(plans[(int)(dp1[0])]);
  }
  else if (strcmp(cmd,"get_x") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for get_x.");

    ARG_GET_PLAN(dp1)
    plhs[0] = mxCreateDoubleMatrix(2, plans[(int)(dp1[0])]->M_total, mxREAL);
    dp2 = mxGetPr(plhs[0]);
    for (i = 0; i < plans[(int)(dp1[0])]->M_total; i++)
    {
      dp2[2*i] = 2*PI*((plans[(int)(dp1[0])]->x[2*i] < 0)?
        (plans[(int)(dp1[0])]->x[2*i] + 1.0):(plans[(int)(dp1[0])]->x[2*i]));
      dp2[2*i+1] = 2*PI*plans[(int)(dp1[0])]->x[2*i+1];
    }
  }
  else if (strcmp(cmd,"get_f") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for get_f.");

    ARG_GET_PLAN(dp1)
    plhs[0] = mxCreateDoubleMatrix(plans[(int)(dp1[0])]->M_total, 1, mxCOMPLEX);
    dp2 = mxGetPr(plhs[0]); dp3 = mxGetPi(plhs[0]);
    for (i = 0; i < plans[(int)(dp1[0])]->M_total; i++)
    {
      dp2[i] = creal(plans[(int)(dp1[0])]->f[i]);
      dp3[i] = cimag(plans[(int)(dp1[0])]->f[i]);
    }
  }
  else if (strcmp(cmd,"get_f_hat") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for get_f_hat.");

    ARG_GET_PLAN(dp1)
    plhs[0] = mxCreateDoubleMatrix(2*plans[(int)(dp1[0])]->N+1,
      plans[(int)(dp1[0])]->N+1, mxCOMPLEX);
    dp2 = mxGetPr(plhs[0]); dp3 = mxGetPi(plhs[0]);

    /*plans[(int)(dp1[0])]->f_hat[NFSFT_INDEX(0, 0, plans[(int)(dp1[0])])] = 1 + 2*_Complex_I;
    plans[(int)(dp1[0])]->f_hat[NFSFT_INDEX(1, -1, plans[(int)(dp1[0])])] = 3 + 4*_Complex_I;
    plans[(int)(dp1[0])]->f_hat[NFSFT_INDEX(1, 0, plans[(int)(dp1[0])])] = 5 + 6*_Complex_I;
    plans[(int)(dp1[0])]->f_hat[NFSFT_INDEX(1, 1, plans[(int)(dp1[0])])] = 7 + 8*_Complex_I;*/

    for (i = 0; i <= plans[(int)(dp1[0])]->N; i++)
    {
      for (j = -i; j <= i; j++)
      {
        dp2[i*(2*plans[(int)(dp1[0])]->N+1)+plans[(int)(dp1[0])]->N+j] =
          creal(plans[(int)(dp1[0])]->f_hat[NFSFT_INDEX(i,j,
            plans[(int)(dp1[0])])]);
        dp3[i*(2*plans[(int)(dp1[0])]->N+1)+plans[(int)(dp1[0])]->N+j] =
          cimag(plans[(int)(dp1[0])]->f_hat[NFSFT_INDEX(i,j,
            plans[(int)(dp1[0])])]);
      }
    }
  }
  else if (strcmp(cmd,"set_x") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 3)
      mexErrMsgTxt("Two additional argument required for set_x.");

    ARG_GET_PLAN(dp1)

    if (mxIsDouble(prhs[2]) != 1 || mxGetNumberOfDimensions(prhs[2]) > 2)
      mexErrMsgTxt("Third argument must be a 2 x M double array");

    if (mxGetM(prhs[2]) != 2 || mxGetN(prhs[2]) != plans[(int)(dp1[0])]->M_total)
      mexErrMsgTxt("Third argument must have correct size.");

    dp2 = mxGetPr(prhs[2]);

    for (i = 0; i < plans[(int)(dp1[0])]->M_total; i++)
    {
      plans[(int)(dp1[0])]->x[2*i] =
        ((dp2[2*i] > PI)?(dp2[2*i] - 2*PI):(dp2[2*i]))/(2*PI);
      plans[(int)(dp1[0])]->x[2*i+1] = dp2[2*i+1]/(2*PI);
    }
  }
  else if (strcmp(cmd,"set_f") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 3)
      mexErrMsgTxt("Two additional argument required for set_f.");

    ARG_GET_PLAN(dp1)

    if (mxIsDouble(prhs[2]) != 1)
      mexErrMsgTxt("Third argument must be a double array");

    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != plans[(int)(dp1[0])]->M_total)
      mexErrMsgTxt("Third argument must have correct size.");

    dp2 = mxGetPr(prhs[2]); dp3 = mxGetPi(prhs[2]);

    for (i = 0; i < plans[(int)(dp1[0])]->M_total; i++)
      plans[(int)(dp1[0])]->f[i] = dp2[i] + ((dp3)?(I*dp3[i]):(0.0));
  }
  else if (strcmp(cmd,"set_f_hat") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 3)
      mexErrMsgTxt("Two additional argument required for set_f_hat.");

    ARG_GET_PLAN(dp1)

    if (mxIsDouble(prhs[2]) != 1)
      mexErrMsgTxt("Third argument must be a double array");

    if (mxGetM(prhs[2]) != 2*plans[(int)(dp1[0])]->N+1 ||
      mxGetN(prhs[2]) != plans[(int)(dp1[0])]->N+1)
      mexErrMsgTxt("Third argument must have correct size.");

    dp2 = mxGetPr(prhs[2]); dp3 = mxGetPi(prhs[2]);

    for (i = 0; i <= plans[(int)(dp1[0])]->N; i++)
    {
      for (j = -i; j <= i; j++)
      {
        plans[(int)(dp1[0])]->f_hat[NFSFT_INDEX(i,j,plans[(int)(dp1[0])])] =
          dp2[i*(2*plans[(int)(dp1[0])]->N+1)+plans[(int)(dp1[0])]->N+j] +
          ((dp3)?(I*dp3[i*(2*plans[(int)(dp1[0])]->N+1)+
                        plans[(int)(dp1[0])]->N+j]):(0.0));
      }
    }
  }
  else if (strcmp(cmd,"display") == 0)
  {
    /* Check number of arguments. */
    if (nrhs != 2)
      mexErrMsgTxt("One additional argument required for display.");

    ARG_GET_PLAN(dp1)
    mexPrintf("Plan %d\n",(int)(dp1[0]));
    mexPrintf("  pointer: %p\n",plans[(int)(dp1[0])]);
    mexPrintf("        N: %d\n",plans[(int)(dp1[0])]->N);
    mexPrintf("  N_total: %d\n",plans[(int)(dp1[0])]->N_total);
    mexPrintf("  M_total: %d\n",plans[(int)(dp1[0])]->M_total);
    mexPrintf("        x: %p\n",plans[(int)(dp1[0])]->x);
    mexPrintf("        f: %p\n",plans[(int)(dp1[0])]->f);
    mexPrintf("    f_hat: %p\n",plans[(int)(dp1[0])]->f_hat);
    mexPrintf("    flags: %d\n",plans[(int)(dp1[0])]->flags);
  }
  else
  {
    mexErrMsgTxt("Unknown command.\n");
  }
}
