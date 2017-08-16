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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "config.h"

#include "nfft3.h"
#include "infft.h"

#define NREPEAT 5

#if defined(_WIN32) || defined(_WIN64)
const char *CMD_CREATEDATASET = "nfsft_benchomp_createdataset.exe";
const char *CMD_DETAIL_SINGLE = "nfsft_benchomp_detail_single.exe";
const char *CMD_DETAIL_THREADS = "nfsft_benchomp_detail_threads.exe";
#else
const char *CMD_CREATEDATASET = "./nfsft_benchomp_createdataset";
const char *CMD_DETAIL_SINGLE = "./nfsft_benchomp_detail_single";
const char *CMD_DETAIL_THREADS = "./nfsft_benchomp_detail_threads";
#endif

static FILE* file_out_tex = NULL;

int get_nthreads_array(int **arr)
{
  int max_threads = X(get_num_threads)();
  int alloc_num = 2;
  int k;
  int ret_number = 0;
  int max_threads_pw2 = (max_threads / 2) * 2 == max_threads ? 1 : 0;

  if (max_threads <= 5)
  {
    *arr = (int*) malloc(max_threads*sizeof(int));
    for (k = 0; k < max_threads; k++)
      *(*arr + k) = k+1;
    return max_threads;
  }

  for (k = 1; k <= max_threads; k*=2, alloc_num++);

  *arr = (int*) malloc(alloc_num*sizeof(int));

  for (k = 1; k <= max_threads; k*=2)
  {
    if (k != max_threads && 2*k > max_threads && max_threads_pw2)
    {
      *(*arr + ret_number) = max_threads/2;
      ret_number++;
    }

    *(*arr + ret_number) = k;
    ret_number++;

    if (k != max_threads && 2*k > max_threads)
    {
      *(*arr + ret_number) = max_threads;
      ret_number++;
      break;
    }
  }

  return ret_number;
} 
  

void check_result_value(const int val, const int ok, const char *msg)
{
  if (val != ok)
  {
    fprintf(stderr, "ERROR %s: %d not %d\n", msg, val, ok);

    exit(1);
  }
}

void run_test_create(int trafo_adjoint, int N, int M)
{
  char cmd[1025];

  snprintf(cmd, 1024, "%s %d %d %d > nfsft_benchomp_test.data", CMD_CREATEDATASET, trafo_adjoint, N, M);
  fprintf(stderr, "%s\n", cmd);
  check_result_value(system(cmd), 0, "createdataset");
}

void run_test_init_output()
{
  FILE *f = fopen("nfsft_benchomp_test.result", "w");
  if (f!= NULL)
    fclose(f);
}

typedef struct
{
  int trafo_adjoint;
  int N;
  int M;
  int m;
  int nfsft_flags;
  int psi_flags;
} s_param;

typedef struct
{
  double avg;
  double min;
  double max;
} s_resval;

typedef struct
{
  int nthreads;
  s_resval resval[6];
} s_result;

typedef struct
{
  s_param param;
  s_result *results;
  int nresults;
} s_testset;

void run_test(s_resval *res, int nrepeat, int m, int nfsft_flags, int psi_flags, int nthreads)
{
  FILE *f;
  char cmd[1025];
  int r,t;
  
  for (t = 0; t < 6; t++)
  {
    res[t].avg = 0.0; res[t].min = 1.0/0.0; res[t].max = 0.0;
  }

  if (nthreads < 2)
    snprintf(cmd, 1024, "%s %d %d %d %d < nfsft_benchomp_test.data > nfsft_benchomp_test.out", CMD_DETAIL_SINGLE, m, nfsft_flags, psi_flags, nrepeat);
  else
    snprintf(cmd, 1024, "%s %d %d %d %d %d < nfsft_benchomp_test.data > nfsft_benchomp_test.out", CMD_DETAIL_THREADS, m, nfsft_flags, psi_flags, nrepeat, nthreads);
  fprintf(stderr, "%s\n", cmd);

  check_result_value(system(cmd), 0, cmd);

  f = fopen("nfsft_benchomp_test.out", "r");
  for (r = 0; r < nrepeat; r++)
  {
    int retval;
    double v[6];
//    FILE *f;
//    check_result_value(system(cmd), 0, cmd);
//    f = fopen("nfsft_benchomp_test.out", "r");
    retval = fscanf(f, "%lg %lg %lg %lg %lg %lg", v, v+1, v+2, v+3, v+4, v+5);
    check_result_value(retval, 6, "read nfsft_benchomp_test.out");
//    fclose(f);
//    fprintf(stderr, "%.3e %.3e %.3e %.3e %.3e %.3e\n", v[0], v[1], v[2], v[3], v[4], v[5]);
    for (t = 0; t < 6; t++)
    {
      res[t].avg += v[t];
      if (res[t].min > v[t])
        res[t].min = v[t];
      if (res[t].max < v[t])
        res[t].max = v[t];
    }
  }
  fclose(f);

  for (t = 0; t < 6; t++)
    res[t].avg /= nrepeat;

  fprintf(stderr, "%d %d: ", nthreads, nrepeat);
  for (t = 0; t < 6; t++)
    fprintf(stderr, "%.3e %.3e %.3e | ", res[t].avg, res[t].min, res[t].max);
  fprintf(stderr, "\n");
}

const char *get_psi_string(int flags)
{
  if (flags & PRE_PSI)
    return "prepsi";
  else if (flags & PRE_ONE_PSI)
    return "unknownPSI";

  return "nopsi";
}
const char *get_sort_string(int flags)
{
  if (flags & NFFT_SORT_NODES)
    return "sorted";

    return "unsorted";
}

const char *get_adjoint_omp_string(int flags)
{
  if (flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    return "blockwise";

    return "";
}

#define MASK_TA (1U<<1)
#define MASK_N (1U<<2)
#define MASK_M (1U<<4)
#define MASK_WINM (1U<<5)
#define MASK_FLAGS_PSI (1U<<6)
#define MASK_FLAGS_SORT (1U<<7)
#define MASK_FLAGS_BW (1U<<8)
#define MASK_FLAGS_FPT (1U<<9)

unsigned int determine_different_parameters(s_testset *testsets, int ntestsets)
{
  int t;
  unsigned int mask = 0;

  if (ntestsets < 2)
    return 0;

  for (t = 1; t < ntestsets; t++)
  {
    if (testsets[t-1].param.trafo_adjoint != testsets[t].param.trafo_adjoint)
      mask |= MASK_TA;
    if (testsets[t-1].param.N != testsets[t].param.N)
      mask |= MASK_N;
    if (testsets[t-1].param.M != testsets[t].param.M)
      mask |= MASK_M;
    if (testsets[t-1].param.m != testsets[t].param.m)
      mask |= MASK_WINM;
    if ((testsets[t-1].param.psi_flags & PRE_ONE_PSI) != (testsets[t].param.psi_flags & PRE_ONE_PSI))
      mask |= MASK_FLAGS_PSI;
    if ((testsets[t-1].param.psi_flags & NFFT_SORT_NODES) != (testsets[t].param.psi_flags & NFFT_SORT_NODES))
      mask |= MASK_FLAGS_SORT;
    if ((testsets[t-1].param.psi_flags & NFFT_OMP_BLOCKWISE_ADJOINT) != (testsets[t].param.psi_flags & NFFT_OMP_BLOCKWISE_ADJOINT))
      mask |= MASK_FLAGS_BW;
    if ((testsets[t-1].param.nfsft_flags & NFSFT_USE_DPT) != (testsets[t].param.nfsft_flags & NFSFT_USE_DPT))
      mask |= MASK_FLAGS_FPT;
  }

  return mask;
}

void get_plot_title(char *outstr, int maxlen, char *hostname, s_param param, unsigned int diff_mask)
{
  unsigned int mask = ~diff_mask;
  int offset = 0;
  int len;

  len = snprintf(outstr, maxlen, "%s", hostname);
  if (len < 0 || len+offset >= maxlen-1) return;
  offset += len;

  if (mask & MASK_TA)
  {
    len = snprintf(outstr+offset, maxlen-offset, " $\\mathrm{NFSFT}%s$", param.trafo_adjoint==0?"":"^\\top");
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_N)
  {
    len = snprintf(outstr+offset, maxlen-offset, " N=%d", param.N);
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_M)
  {
    len = snprintf(outstr+offset, maxlen-offset, " M=%d", param.M);
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_WINM)
  {
    len = snprintf(outstr+offset, maxlen-offset, " m=%d", param.m);
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_FLAGS_PSI)
  {
    len = snprintf(outstr+offset, maxlen-offset, " %s", get_psi_string(param.psi_flags));
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_FLAGS_SORT)
  {
    len = snprintf(outstr+offset, maxlen-offset, " %s", get_sort_string(param.psi_flags));
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if ((mask & MASK_FLAGS_BW) && strlen(get_adjoint_omp_string(param.psi_flags)) > 0)
  {
    len = snprintf(outstr+offset, maxlen-offset, " %s", get_adjoint_omp_string(param.psi_flags));
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_FLAGS_FPT)
  {
    len = snprintf(outstr+offset, maxlen-offset, param.nfsft_flags & NFSFT_USE_DPT ? " DPT" : "");
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

}

void print_output_speedup_total_tref(FILE *out, s_testset *testsets, int ntestsets, int use_tref, double tref)
{
  int i, t;
  char hostname[1025];
  char plottitle[1025];
  unsigned int diff_mask = determine_different_parameters(testsets, ntestsets);

#ifdef HAVE_GETHOSTNAME
  if (gethostname(hostname, 1024) != 0)
#endif
    strncpy(hostname, "unnamed", 1024);

  get_plot_title(plottitle, 1024, hostname, testsets[0].param, diff_mask);

  fprintf(out, "\\begin{tikzpicture}\n");
  fprintf(out, "\\begin{axis}[");
  fprintf(out, "width=0.9\\textwidth, height=0.6\\textwidth, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Speedup, xtick=data, legend style={ legend pos = north west, legend columns=1}, ymajorgrids=true, yminorgrids=true, minor y tick num=4, ");
  fprintf(out, " title={%s}", plottitle);
  fprintf(out, " ]\n");

  for (t = 0; t < ntestsets; t++)
  {
    s_testset testset = testsets[t];
    fprintf(stderr, "%s $\\mathrm{NFSFT}%s$ N=%d M=%d m=%d %s %s %s}", hostname, testset.param.trafo_adjoint==0?"":"^\\top", testset.param.N, testset.param.M, testset.param.m, get_psi_string(testset.param.psi_flags), get_sort_string(testset.param.psi_flags), get_adjoint_omp_string(testset.param.psi_flags));
    fprintf(stderr, "\n");

    fprintf(out, "\\addplot coordinates {");
    for (i = 0; i < testset.nresults; i++)
      if (use_tref == 1)
        fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, tref/testset.results[i].resval[5].avg);
      else
        fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, testset.results[0].resval[5].avg/testset.results[i].resval[5].avg);
    fprintf(out, "};\n");

    for (i = 0; i < testset.nresults; i++)
      if (use_tref == 1)
        fprintf(stderr, "%d:%.3f  ", testset.results[i].nthreads, tref/testset.results[i].resval[5].avg);
      else
        fprintf(stderr, "%d:%.3f  ", testset.results[i].nthreads, testset.results[0].resval[5].avg/testset.results[i].resval[5].avg);
    fprintf(stderr, "\n\n");
  }

  fprintf(out, "\\legend{{");
  for (t = 0; t < ntestsets; t++)
  {
    char title[256];
    if (t > 0)
      fprintf(out, "},{");
    get_plot_title(title, 255, "", testsets[t].param, ~(diff_mask));
    fprintf(out, "%s", title);
  }
  fprintf(out, "}}\n");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void print_output_speedup_total(FILE *out, s_testset *testsets, int ntestsets, int use_tref)
{
  double tref = 1.0/0.0;
  int t, k;

  if (use_tref == 1)
    for (t = 0; t < ntestsets; t++)
      for (k = 0; k < testsets[t].nresults; k++)
        if (testsets[t].results[k].nthreads == 1 && testsets[t].results[k].resval[5].avg < tref)
          tref = testsets[t].results[k].resval[5].avg;

  print_output_speedup_total_tref(out, testsets, ntestsets, use_tref, tref);
}

void print_output_histo_PENRT(FILE *out, s_testset testset)
{
  int i, size = testset.nresults;
  char hostname[1025];

#ifdef HAVE_GETHOSTNAME
  if (gethostname(hostname, 1024) != 0)
#endif
    strncpy(hostname, "unnamed", 1024);

  fprintf(out, "\\begin{tikzpicture}\n");
  fprintf(out, "\\begin{axis}[");
  fprintf(out, "width=0.9\\textwidth, height=0.6\\textwidth, ");
  fprintf(out, "symbolic x coords={");
  for (i = 0; i < size; i++)
    if (i > 0)
      fprintf(out, ",%d", testset.results[i].nthreads);
    else
      fprintf(out, "%d", testset.results[i].nthreads);

  fprintf(out, "}, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Time in s, xtick=data, legend style={legend columns=-1}, ybar, bar width=7pt, ymajorgrids=true, yminorgrids=true, minor y tick num=1, ");
  fprintf(out, " title={%s $\\mathrm{NFSFT}%s$ N=%d M=%d m=%d %s %s %s}", hostname, testset.param.trafo_adjoint==0?"":"^\\top", testset.param.N, testset.param.M, testset.param.m, get_psi_string(testset.param.psi_flags), get_sort_string(testset.param.psi_flags), get_adjoint_omp_string(testset.param.psi_flags));
  fprintf(out, " ]\n");
  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, testset.results[i].resval[1].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, testset.results[i].resval[2].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, testset.results[i].resval[3].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, testset.results[i].resval[0].avg + testset.results[i].resval[4].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, testset.results[i].resval[5].avg);
  fprintf(out, "};\n");
  fprintf(out, "\\legend{%s,%s,$\\mathrm{NFFT}%s$,rest,total}\n", testset.param.nfsft_flags & NFSFT_USE_DPT ? "DPT" : "FPT", testset.param.trafo_adjoint==0?"c2e":"$\\mathrm{c2e}^\\top$", testset.param.trafo_adjoint==0?"":"^\\top");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void run_testset(s_testset *testset, int trafo_adjoint, int N, int M, int m, int nfsft_flags, int psi_flags, int *nthreads_array, int n_threads_array_size)
{
  int i;
  testset->param.trafo_adjoint = trafo_adjoint;
  testset->param.N = N;
  testset->param.M = M;
  testset->param.m = m;
  testset->param.nfsft_flags = nfsft_flags;
  testset->param.psi_flags = psi_flags;

  testset->results = (s_result*) malloc(n_threads_array_size*sizeof(s_result));
  testset->nresults = n_threads_array_size;

  run_test_create(testset->param.trafo_adjoint, testset->param.N, testset->param.M);
  for (i = 0; i < n_threads_array_size; i++)
  {
    testset->results[i].nthreads = nthreads_array[i];
    run_test(testset->results[i].resval, NREPEAT, testset->param.m, testset->param.nfsft_flags, testset->param.psi_flags, testset->results[i].nthreads = nthreads_array[i]);
  }

}

void test1(int *nthreads_array, int n_threads_array_size, int m)
{
  s_testset testsets[4];

  run_testset(&testsets[0], 0, 1024, 1000000, m, 0, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_PENRT(file_out_tex, testsets[0]);
#endif

  run_testset(&testsets[1], 1, 1024, 1000000, m, 0, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_PENRT(file_out_tex, testsets[1]);
#endif

  print_output_speedup_total(file_out_tex, testsets, 2, 0);

  run_testset(&testsets[2], 0, 1024, 1000000, m, NFSFT_USE_DPT, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_PENRT(file_out_tex, testsets[2]);
#endif

  run_testset(&testsets[3], 1, 1024, 1000000, m, NFSFT_USE_DPT, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_PENRT(file_out_tex, testsets[3]);
#endif

  print_output_speedup_total(file_out_tex, testsets+2, 2, 0);
}

int main(int argc, char** argv)
{
  int *nthreads_array;
  int n_threads_array_size = get_nthreads_array(&nthreads_array);
  int k;

#if !(defined MEASURE_TIME && defined MEASURE_TIME_FFTW)
  fprintf(stderr, "WARNING: Detailed time measurements for NFSFT are not activated.\n");
  fprintf(stderr, "For more detailed plots, please re-run the configure script with options\n");
  fprintf(stderr, "--enable-measure-time --enable-measure-time-fftw --enable-nfsft --enable-openmp\n");
  fprintf(stderr, "and run \"make clean all\"\n\n");
#endif

  for (k = 0; k < n_threads_array_size; k++)
    fprintf(stderr, "%d ", nthreads_array[k]);
  fprintf(stderr, "\n");

  file_out_tex = fopen("nfsft_benchomp_results_plots.tex", "w");

  test1(nthreads_array, n_threads_array_size, 2);
  test1(nthreads_array, n_threads_array_size, 4);
  test1(nthreads_array, n_threads_array_size, 6);
  test1(nthreads_array, n_threads_array_size, 8);

  fclose(file_out_tex);

  return 0;
}
