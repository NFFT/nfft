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
const char *CMD_CREATEDATASET = "nfft_benchomp_createdataset.exe";
const char *CMD_DETAIL_SINGLE = "nfft_benchomp_detail_single.exe";
const char *CMD_DETAIL_THREADS = "nfft_benchomp_detail_threads.exe";
#else
const char *CMD_CREATEDATASET = "./nfft_benchomp_createdataset";
const char *CMD_DETAIL_SINGLE = "./nfft_benchomp_detail_single";
const char *CMD_DETAIL_THREADS = "./nfft_benchomp_detail_threads";
#endif

static FILE* file_out_tex = NULL;

int get_nthreads_array(int **arr)
{
  int max_threads = NFFT(get_num_threads)();
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

void run_test_create(int d, int trafo_adjoint, int N, int M, double sigma)
{
  char cmd[1025];

  if (d==1)
    snprintf(cmd, 1024, "%s %d %d %d %d %lg > nfft_benchomp_test.data", CMD_CREATEDATASET, d, trafo_adjoint, N, M, sigma);
  else if (d==2)  
    snprintf(cmd, 1024, "%s %d %d %d %d %d %lg > nfft_benchomp_test.data", CMD_CREATEDATASET, d, trafo_adjoint, N, N, M, sigma);
  else if (d==3)  
    snprintf(cmd, 1024, "%s %d %d %d %d %d %d %lg > nfft_benchomp_test.data", CMD_CREATEDATASET, d, trafo_adjoint, N, N, N, M, sigma);
  else if (d==4)  
    snprintf(cmd, 1024, "%s %d %d %d %d %d %d %d %lg > nfft_benchomp_test.data", CMD_CREATEDATASET, d, trafo_adjoint, N, N, N, N, M, sigma);
  else
    exit(1);
  fprintf(stderr, "%s\n", cmd);
  check_result_value(system(cmd), 0, "createdataset");
}

void run_test_init_output()
{
  FILE *f = fopen("nfft_benchomp_test.result", "w");
  if (f!= NULL)
    fclose(f);
}

typedef struct
{
  int d;
  int trafo_adjoint;
  int N;
  int M;
  double sigma;
  int m;
  int flags;
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

void run_test(s_resval *res, int nrepeat, int m, int flags, int nthreads)
{
  char cmd[1025];
  int r,t;
  
  for (t = 0; t < 6; t++)
  {
    res[t].avg = 0.0; res[t].min = 1.0/0.0; res[t].max = 0.0;
  }

  if (nthreads < 2)
    snprintf(cmd, 1024, "%s %d %d < nfft_benchomp_test.data > nfft_benchomp_test.out", CMD_DETAIL_SINGLE, m, flags);
  else
    snprintf(cmd, 1024, "%s %d %d %d < nfft_benchomp_test.data > nfft_benchomp_test.out", CMD_DETAIL_THREADS, m, flags, nthreads);
  fprintf(stderr, "%s\n", cmd);
  check_result_value(system(cmd), 0, cmd);

  for (r = 0; r < nrepeat; r++)
  {
    int retval;
    double v[6];
    FILE *f;
    check_result_value(system(cmd), 0, cmd);
    f = fopen("nfft_benchomp_test.out", "r");
    retval = fscanf(f, "%lg %lg %lg %lg %lg %lg", v, v+1, v+2, v+3, v+4, v+5);
    check_result_value(retval, 6, "read nfft_benchomp_test.out");
    fclose(f);

    for (t = 0; t < 6; t++)
    {
      res[t].avg += v[t];
      if (res[t].min > v[t])
        res[t].min = v[t];
      if (res[t].max < v[t])
        res[t].max = v[t];
    }
  }

  for (t = 0; t < 6; t++)
    res[t].avg /= nrepeat;

  fprintf(stderr, "%d %d: ", nthreads, nrepeat);
  for (t = 0; t < 6; t++)
    fprintf(stderr, "%.3e %.3e %.3e | ", res[t].avg, res[t].min, res[t].max);
  fprintf(stderr, "\n");
}

const char *get_psi_string(int flags)
{
  if (flags & PRE_ONE_PSI)
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

#define MASK_D (1U<<0)
#define MASK_TA (1U<<1)
#define MASK_N (1U<<2)
#define MASK_SIGMA (1U<<3)
#define MASK_M (1U<<4)
#define MASK_WINM (1U<<5)
#define MASK_FLAGS_PSI (1U<<6)
#define MASK_FLAGS_SORT (1U<<7)
#define MASK_FLAGS_BW (1U<<8)

unsigned int determine_different_parameters(s_testset *testsets, int ntestsets)
{
  int t;
  unsigned int mask = 0;

  if (ntestsets < 2)
    return 0;

  for (t = 1; t < ntestsets; t++)
  {
    if (testsets[t-1].param.d != testsets[t].param.d)
      mask |= MASK_D;
    if (testsets[t-1].param.trafo_adjoint != testsets[t].param.trafo_adjoint)
      mask |= MASK_TA;
    if (testsets[t-1].param.N != testsets[t].param.N)
      mask |= MASK_N;
    if (testsets[t-1].param.sigma != testsets[t].param.sigma)
      mask |= MASK_SIGMA;
    if (testsets[t-1].param.M != testsets[t].param.M)
      mask |= MASK_M;
    if (testsets[t-1].param.m != testsets[t].param.m)
      mask |= MASK_WINM;
    if ((testsets[t-1].param.flags & PRE_ONE_PSI) != (testsets[t].param.flags & PRE_ONE_PSI))
      mask |= MASK_FLAGS_PSI;
    if ((testsets[t-1].param.flags & NFFT_SORT_NODES) != (testsets[t].param.flags & NFFT_SORT_NODES))
      mask |= MASK_FLAGS_SORT;
    if ((testsets[t-1].param.flags & NFFT_OMP_BLOCKWISE_ADJOINT) != (testsets[t].param.flags & NFFT_OMP_BLOCKWISE_ADJOINT))
      mask |= MASK_FLAGS_BW;
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

  if (mask & MASK_D)
  {
    len = snprintf(outstr+offset, maxlen-offset, " %dd", param.d);
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_TA)
  {
    len = snprintf(outstr+offset, maxlen-offset, " $\\mathrm{NFFT}%s$", param.trafo_adjoint==0?"":"^\\top");
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_N)
  {
    len = snprintf(outstr+offset, maxlen-offset, " N=%d", param.N);
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_SIGMA)
  {
    len = snprintf(outstr+offset, maxlen-offset, " N=%g", param.sigma);
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
    len = snprintf(outstr+offset, maxlen-offset, " %s", get_psi_string(param.flags));
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if (mask & MASK_FLAGS_SORT)
  {
    len = snprintf(outstr+offset, maxlen-offset, " %s", get_sort_string(param.flags));
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }

  if ((mask & MASK_FLAGS_BW) && strlen(get_adjoint_omp_string(param.flags)) > 0)
  {
    len = snprintf(outstr+offset, maxlen-offset, " %s", get_adjoint_omp_string(param.flags));
    if (len < 0 || len+offset >= maxlen-1) return;
    offset += len;
  }
}

void print_output_speedup_total_tref(FILE *out, s_testset *testsets, int ntestsets, double tref)
{
  int i, t;
  char hostname[1025];
  char plottitle[1025];
  unsigned int diff_mask = determine_different_parameters(testsets, ntestsets);

#ifdef HAVE_GETHOSTNAME
  if (gethostname(hostname, 1024) != 0)
#endif
    strncpy(hostname, "unnamed", 1024);

  get_plot_title(plottitle, 1024, hostname, testsets[0].param, diff_mask | MASK_FLAGS_SORT);

  fprintf(out, "\\begin{tikzpicture}\n");
  fprintf(out, "\\begin{axis}[");
  fprintf(out, "width=0.9\\textwidth, height=0.6\\textwidth, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Speedup, xtick=data, legend style={ legend pos = north west, legend columns=1}, ymajorgrids=true, yminorgrids=true, minor y tick num=4, ");
  fprintf(out, " title={%s}", plottitle);
  fprintf(out, " ]\n");

  for (t = 0; t < ntestsets; t++)
  {
    s_testset testset = testsets[t];
    fprintf(stderr, "%s %dd $\\mathrm{NFFT}%s$ N=%d $\\sigma$=%g M=%d m=%d %s %s %s}", hostname, testset.param.d, testset.param.trafo_adjoint==0?"":"^\\top", testset.param.N, testset.param.sigma, testset.param.M, testset.param.m, get_psi_string(testset.param.flags), get_sort_string(testset.param.flags), get_adjoint_omp_string(testset.param.flags));
    fprintf(stderr, "\n");

    fprintf(out, "\\addplot coordinates {");
    for (i = 0; i < testset.nresults; i++)
      fprintf(out, "(%d, %.6e) ", testset.results[i].nthreads, tref/testset.results[i].resval[5].avg);
    fprintf(out, "};\n");

    for (i = 0; i < testset.nresults; i++)
    {
      fprintf(stderr, "%d:%.3f  ", testset.results[i].nthreads, tref/testset.results[i].resval[5].avg);
    }
    fprintf(stderr, "\n\n");
  }

  fprintf(out, "\\legend{{");
  for (t = 0; t < ntestsets; t++)
  {
    char title[256];
    if (t > 0)
      fprintf(out, "},{");
    get_plot_title(title, 255, "", testsets[t].param, ~(diff_mask | MASK_FLAGS_SORT));
    fprintf(out, "%s", title);
  }
  fprintf(out, "}}\n");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void print_output_speedup_total(FILE *out, s_testset *testsets, int ntestsets)
{
  double tref = 1.0/0.0;
  int t, k;

  for (t = 0; t < ntestsets; t++)
    for (k = 0; k < testsets[t].nresults; k++)
      if (testsets[t].results[k].nthreads == 1 && testsets[t].results[k].resval[5].avg < tref)
        tref = testsets[t].results[k].resval[5].avg;

  print_output_speedup_total_tref(out, testsets, ntestsets, tref);
}

void print_output_histo_DFBRT(FILE *out, s_testset testset)
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
fprintf(stderr, "FLAGS: %d\n", testset.param.flags);

  fprintf(out, "}, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Time in s, xtick=data, legend style={legend columns=-1}, ybar, bar width=7pt, ymajorgrids=true, yminorgrids=true, minor y tick num=1, ");
  fprintf(out, " title={%s %dd $\\mathrm{NFFT}%s$ N=%d $\\sigma$=%g M=%d m=%d %s %s %s}", hostname, testset.param.d, testset.param.trafo_adjoint==0?"":"^\\top", testset.param.N, testset.param.sigma, testset.param.M, testset.param.m, get_psi_string(testset.param.flags), get_sort_string(testset.param.flags), get_adjoint_omp_string(testset.param.flags));
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
  fprintf(out, "\\legend{D,F,B,rest,total}\n");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void run_testset(s_testset *testset, int d, int trafo_adjoint, int N, int M, double sigma, int m, int flags, int *nthreads_array, int n_threads_array_size)
{
  int i;
  testset->param.d = d;
  testset->param.trafo_adjoint = trafo_adjoint;
  testset->param.N = N;
  testset->param.M = M;
  testset->param.sigma = sigma;
  testset->param.m = m;
  testset->param.flags = flags;

  testset->results = (s_result*) malloc(n_threads_array_size*sizeof(s_result));
  testset->nresults = n_threads_array_size;

  run_test_create(testset->param.d, testset->param.trafo_adjoint, testset->param.N, testset->param.M, testset->param.sigma);
  for (i = 0; i < n_threads_array_size; i++)
  {
    testset->results[i].nthreads = nthreads_array[i];
    run_test(testset->results[i].resval, NREPEAT, testset->param.m, testset->param.flags, testset->results[i].nthreads = nthreads_array[i]);
  }

}

void test1(int *nthreads_array, int n_threads_array_size, int m)
{
  s_testset testsets[15];

  run_testset(&testsets[0], 1, 0, 2097152, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[0]);
#endif

  run_testset(&testsets[1], 1, 0, 2097152, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[1]);
#endif

  print_output_speedup_total(file_out_tex, testsets, 2);

  run_testset(&testsets[2], 1, 1, 2097152, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[2]);
#endif

  run_testset(&testsets[3], 1, 1, 2097152, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[3]);
#endif

  run_testset(&testsets[4], 1, 1, 2097152, 2097152, 2.0, m, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[4]);
#endif

  print_output_speedup_total(file_out_tex, testsets+2, 3);

  run_testset(&testsets[5], 2, 0, 1024, 1048576, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[5]);
#endif

  run_testset(&testsets[6], 2, 0, 1024, 1048576, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[6]);
#endif

  print_output_speedup_total(file_out_tex, testsets+5, 2);

  run_testset(&testsets[7], 2, 1, 1024, 1048576, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[7]);
#endif

  run_testset(&testsets[8], 2, 1, 1024, 1048576, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[8]);
#endif

  run_testset(&testsets[9], 2, 1, 1024, 1048576, 2.0, m, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[9]);
#endif

  print_output_speedup_total(file_out_tex, testsets+7, 3);

  run_testset(&testsets[10], 3, 0, 128, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[10]);
#endif

  run_testset(&testsets[11], 3, 0, 128, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[11]);
#endif

  print_output_speedup_total(file_out_tex, testsets+10, 2);

  run_testset(&testsets[12], 3, 1, 128, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[12]);
#endif

  run_testset(&testsets[13], 3, 1, 128, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[13]);
#endif

  run_testset(&testsets[14], 3, 1, 128, 2097152, 2.0, m, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[14]);
#endif

  print_output_speedup_total(file_out_tex, testsets+12, 3);

}

void test2(int *nthreads_array, int n_threads_array_size, int m)
{
  s_testset testsets[15];

  run_testset(&testsets[0], 1, 0, 16777216, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[0]);
#endif

  run_testset(&testsets[1], 1, 0, 16777216, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[1]);
#endif

  print_output_speedup_total(file_out_tex, testsets, 2);

  run_testset(&testsets[2], 1, 1, 16777216, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[2]);
#endif

  run_testset(&testsets[3], 1, 1, 16777216, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[3]);
#endif

  run_testset(&testsets[4], 1, 1, 16777216, 2097152, 2.0, m, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[4]);
#endif

  print_output_speedup_total(file_out_tex, testsets+2, 3);

  run_testset(&testsets[5], 2, 0, 4096, 1048576, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[5]);
#endif

  run_testset(&testsets[6], 2, 0, 4096, 1048576, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[6]);
#endif

  print_output_speedup_total(file_out_tex, testsets+5, 2);

  run_testset(&testsets[7], 2, 1, 4096, 1048576, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[7]);
#endif

  run_testset(&testsets[8], 2, 1, 4096, 1048576, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[8]);
#endif

  run_testset(&testsets[9], 2, 1, 4096, 1048576, 2.0, m, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[9]);
#endif

  print_output_speedup_total(file_out_tex, testsets+7, 3);

  run_testset(&testsets[10], 3, 0, 256, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[10]);
#endif

  run_testset(&testsets[11], 3, 0, 256, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[11]);
#endif

  print_output_speedup_total(file_out_tex, testsets+10, 2);

  run_testset(&testsets[12], 3, 1, 256, 2097152, 2.0, m, 0, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[12]);
#endif

  run_testset(&testsets[13], 3, 1, 256, 2097152, 2.0, m, NFFT_SORT_NODES, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[13]);
#endif

  run_testset(&testsets[14], 3, 1, 256, 2097152, 2.0, m, NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT, nthreads_array, n_threads_array_size);
#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  print_output_histo_DFBRT(file_out_tex, testsets[14]);
#endif

  print_output_speedup_total(file_out_tex, testsets+12, 3);

}

int main(int argc, char** argv)
{
  int *nthreads_array;
  int n_threads_array_size = get_nthreads_array(&nthreads_array);
  int k;

#if !(defined MEASURE_TIME && defined MEASURE_TIME_FFTW)
  fprintf(stderr, "WARNING: Detailed time measurements for NFFT are not activated.\n");
  fprintf(stderr, "For more detailed plots, please re-run the configure script with options\n");
  fprintf(stderr, "--enable-measure-time --enable-measure-time-fftw --enable-openmp\n");
  fprintf(stderr, "and run \"make clean all\"\n\n");
#endif

  for (k = 0; k < n_threads_array_size; k++)
    fprintf(stderr, "%d ", nthreads_array[k]);
  fprintf(stderr, "\n");

  file_out_tex = fopen("nfft_benchomp_results_plots.tex", "w");

  test1(nthreads_array, n_threads_array_size, 2);
  test1(nthreads_array, n_threads_array_size, 4);
  test1(nthreads_array, n_threads_array_size, 6);
//  test2(nthreads_array, n_threads_array_size, 2);

  fclose(file_out_tex);

  return 0;
}
