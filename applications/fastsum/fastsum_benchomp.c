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
const char *CMD_CREATEDATASET = "fastsum_benchomp_createdataset.exe";
const char *CMD_DETAIL_SINGLE = "fastsum_benchomp_detail_single.exe";
const char *CMD_DETAIL_THREADS = "fastsum_benchomp_detail_threads.exe";
#else
const char *CMD_CREATEDATASET = "./fastsum_benchomp_createdataset";
const char *CMD_DETAIL_SINGLE = "./fastsum_benchomp_detail_single";
const char *CMD_DETAIL_THREADS = "./fastsum_benchomp_detail_threads";
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
    *arr = (int*) NFFT(malloc)((size_t) (max_threads) * sizeof(int));
    for (k = 0; k < max_threads; k++)
      *(*arr + k) = k + 1;
    return max_threads;
  }

  for (k = 1; k <= max_threads; k *= 2, alloc_num++)
    ;

  *arr = (int*) NFFT(malloc)((size_t)(alloc_num) * sizeof(int));

  for (k = 1; k <= max_threads; k *= 2)
  {
    if (k != max_threads && 2 * k > max_threads && max_threads_pw2)
    {
      *(*arr + ret_number) = max_threads / 2;
      ret_number++;
    }

    *(*arr + ret_number) = k;
    ret_number++;

    if (k != max_threads && 2 * k > max_threads)
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

    exit(EXIT_FAILURE);
  }
}

void run_test_create(int d, int L, int M)
{
  char cmd[1025];

  snprintf(cmd, 1024,
      "%s %d %d %d > fastsum_benchomp_test.data",
      CMD_CREATEDATASET, d, L, M);
  fprintf(stderr, "%s\n", cmd);
  check_result_value(system(cmd), 0, "createdataset");
}

void run_test_init_output()
{
  FILE *f = fopen("fastsum_benchomp_test.result", "w");
  if (f != NULL)
    fclose(f);
}

typedef struct
{
  int d;
  int L;
  int M;
  int n;
  int m;
  int p;
  char *kernel_name;
  R c;
  R eps_I;
  R eps_B;
} s_param;

typedef struct
{
  R avg;
  R min;
  R max;
} s_resval;

typedef struct
{
  int nthreads;
  s_resval resval[16];
} s_result;

typedef struct
{
  s_param param;
  s_result *results;
  int nresults;
} s_testset;

void run_test(s_resval *res, int nrepeat, int n, int m, int p,
    char *kernel_name, R c, R eps_I, R eps_B, int nthreads)
{
  char cmd[1025];
  int r, t;

  for (t = 0; t < 16; t++)
  {
    res[t].avg = K(0.0);
    res[t].min = K(1.0) / K(0.0);
    res[t].max = K(0.0);
  }

  if (nthreads < 2)
    snprintf(cmd, 1024,
        "%s %d %d %d %s " __FR__ " " __FR__ " " __FR__ " < fastsum_benchomp_test.data > fastsum_benchomp_test.out",
        CMD_DETAIL_SINGLE, n, m, p, kernel_name, c, eps_I, eps_B);
  else
    snprintf(cmd, 1024,
        "%s %d %d %d %s " __FR__ " " __FR__ " " __FR__ " %d < fastsum_benchomp_test.data > fastsum_benchomp_test.out",
        CMD_DETAIL_THREADS, n, m, p, kernel_name, c, eps_I, eps_B, nthreads);
  fprintf(stderr, "%s\n", cmd);
  check_result_value(system(cmd), 0, cmd);

  for (r = 0; r < nrepeat; r++)
  {
    int retval;
    R v[16];
    FILE *f;
    check_result_value(system(cmd), 0, cmd);
    f = fopen("fastsum_benchomp_test.out", "r");
    retval = fscanf(f,
        "" __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ " " __FR__ "", v,
        v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7, v + 8, v + 9, v + 10,
        v + 11, v + 12, v + 13, v + 14, v + 15);
    check_result_value(retval, 16, "read fastsum_benchomp_test.out");
    fclose(f);

    for (t = 0; t < 16; t++)
    {
      res[t].avg += v[t];
      if (res[t].min > v[t])
        res[t].min = v[t];
      if (res[t].max < v[t])
        res[t].max = v[t];
    }
  }

  for (t = 0; t < 16; t++)
    res[t].avg /= (R)(nrepeat);

  fprintf(stderr, "%d %d: ", nthreads, nrepeat);
  for (t = 0; t < 16; t++)
    fprintf(stderr, "%.3" __FES__ " %.3" __FES__ " %.3" __FES__ " | ", res[t].avg, res[t].min, res[t].max);
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
  if (flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    return "";

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

#define MASK_FSUM_D (1U<<0)
#define MASK_FSUM_L (1U<<1)
#define MASK_FSUM_M (1U<<2)
#define MASK_FSUM_MULTIBW (1U<<3)
#define MASK_FSUM_WINM (1U<<4)
#define MASK_FSUM_P (1U<<5)
#define MASK_FSUM_KERNEL (1U<<6)
#define MASK_FSUM_EPSI (1U<<7)
#define MASK_FSUM_EPSB (1U<<8)

unsigned int fastsum_determine_different_parameters(s_testset *testsets,
    int ntestsets)
{
  int t;
  unsigned int mask = 0;

  if (ntestsets < 2)
    return 0;

  for (t = 1; t < ntestsets; t++)
  {
    if (testsets[t - 1].param.d != testsets[t].param.d)
      mask |= MASK_FSUM_D;
    if (testsets[t - 1].param.L != testsets[t].param.L)
      mask |= MASK_FSUM_L;
    if (testsets[t - 1].param.M != testsets[t].param.M)
      mask |= MASK_FSUM_M;
    if (testsets[t - 1].param.n != testsets[t].param.n)
      mask |= MASK_FSUM_MULTIBW;
    if (testsets[t - 1].param.m != testsets[t].param.m)
      mask |= MASK_FSUM_WINM;
    if (testsets[t - 1].param.p != testsets[t].param.p)
      mask |= MASK_FSUM_P;
    if (strcmp(testsets[t - 1].param.kernel_name, testsets[t].param.kernel_name)
        != 0)
      mask |= MASK_FSUM_KERNEL;
    if (testsets[t - 1].param.eps_I != testsets[t].param.eps_I)
      mask |= MASK_FSUM_EPSI;
    if (testsets[t - 1].param.eps_B != testsets[t].param.eps_B)
      mask |= MASK_FSUM_EPSB;
  }

  return mask;
}

void strEscapeUnderscore(char *dst, char *src, int maxlen)
{
  int i = 0;
  int len;
  int offset = 0;

  while (src[i] != '\0' && len + offset < maxlen - 1)
  {
    if (src[i] == '_')
      len = snprintf(dst + offset, maxlen - offset, "\\_{}");
    else
      len = snprintf(dst + offset, maxlen - offset, "%c", src[i]);
    offset += len;
    i++;
  }
}

void fastsum_get_plot_title_minus_indep(char *outstr, int maxlen,
    char *hostname, s_param param, unsigned int diff_mask)
{
  unsigned int mask = ~diff_mask;
  int offset = 0;
  int len;

  len = snprintf(outstr, maxlen, "%s", hostname);
  if (len < 0 || len + offset >= maxlen - 1)
    return;
  offset += len;

  if (mask & MASK_FSUM_D)
  {
    len = snprintf(outstr + offset, maxlen - offset, " %dd fastsum", param.d);
    if (len < 0 || len + offset >= maxlen - 1)
      return;
    offset += len;
  }

  if ((mask & (MASK_FSUM_L | MASK_FSUM_M)) && param.L == param.M)
  {
    len = snprintf(outstr + offset, maxlen - offset, " L=M=%d", param.L);
    if (len < 0 || len + offset >= maxlen - 1)
      return;
    offset += len;
  }
  else
  {
    if (mask & MASK_FSUM_L)
    {
      len = snprintf(outstr + offset, maxlen - offset, " L=%d", param.L);
      if (len < 0 || len + offset >= maxlen - 1)
        return;
      offset += len;
    }

    if (mask & MASK_FSUM_M)
    {
      len = snprintf(outstr + offset, maxlen - offset, " M=%d", param.M);
      if (len < 0 || len + offset >= maxlen - 1)
        return;
      offset += len;
    }
  }

  if (mask & MASK_FSUM_MULTIBW)
  {
    len = snprintf(outstr + offset, maxlen - offset, " n=%d", param.n);
    if (len < 0 || len + offset >= maxlen - 1)
      return;
    offset += len;
  }

  if (mask & MASK_FSUM_WINM)
  {
    len = snprintf(outstr + offset, maxlen - offset, " m=%d", param.m);
    if (len < 0 || len + offset >= maxlen - 1)
      return;
    offset += len;
  }

  if (mask & MASK_FSUM_P)
  {
    len = snprintf(outstr + offset, maxlen - offset, " p=%d", param.p);
    if (len < 0 || len + offset >= maxlen - 1)
      return;
    offset += len;
  }

  if (mask & MASK_FSUM_KERNEL)
  {
    char tmp[maxlen];
    strEscapeUnderscore(tmp, param.kernel_name, maxlen);

    len = snprintf(outstr + offset, maxlen - offset, " %s", tmp);
    if (len < 0 || len + offset >= maxlen - 1)
      return;
    offset += len;
  }

  if ((mask & (MASK_FSUM_EPSI | MASK_FSUM_EPSB)) && param.eps_I == param.eps_B)
  {
    len = snprintf(outstr + offset, maxlen - offset,
        " $\\varepsilon_\\mathrm{I}$=$\\varepsilon_\\mathrm{B}$=%" __FGS__ "",
        param.eps_I);
    if (len < 0 || len + offset >= maxlen - 1)
      return;
    offset += len;
  }
  else
  {
    if (mask & MASK_FSUM_EPSI)
    {
      len = snprintf(outstr + offset, maxlen - offset,
          " $\\varepsilon_\\mathrm{I}$=%" __FGS__ "", param.eps_I);
      if (len < 0 || len + offset >= maxlen - 1)
        return;
      offset += len;
    }

    if (mask & MASK_FSUM_EPSB)
    {
      len = snprintf(outstr + offset, maxlen - offset,
          " $\\varepsilon_\\mathrm{B}$=%" __FGS__ "", param.eps_B);
      if (len < 0 || len + offset >= maxlen - 1)
        return;
      offset += len;
    }
  }
}

void nfft_adjoint_print_output_histo_DFBRT(FILE *out, s_testset testset)
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

  fprintf(out,
      "}, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Time in s, xtick=data, legend style={legend columns=-1}, ybar, bar width=7pt, ymajorgrids=true, yminorgrids=true, minor y tick num=1, ");
  fprintf(out,
      " title={%s %dd $\\textrm{NFFT}^\\top$ N=%d $\\sigma$=2 M=%d m=%d prepsi sorted}",
      hostname, testset.param.d, testset.param.n, testset.param.M,
      testset.param.m);
  fprintf(out, " ]\n");
  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[10].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[11].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[12].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[1].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[4].avg + testset.results[i].resval[1].avg);
  fprintf(out, "};\n");
  fprintf(out,
      "\\legend{D,$\\textrm{F}^\\top$,$\\textrm{B}^\\top$,prepsi,total}\n");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void nfft_trafo_print_output_histo_DFBRT(FILE *out, s_testset testset)
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

  fprintf(out,
      "}, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Time in s, xtick=data, legend style={legend columns=-1}, ybar, bar width=7pt, ymajorgrids=true, yminorgrids=true, minor y tick num=1, ");
  fprintf(out,
      " title={%s %dd $\\textrm{NFFT}$ N=%d $\\sigma$=2 M=%d m=%d prepsi sorted}",
      hostname, testset.param.d, testset.param.n, testset.param.M,
      testset.param.m);
  fprintf(out, " ]\n");
  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[13].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[14].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[15].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[2].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[6].avg + testset.results[i].resval[2].avg);
  fprintf(out, "};\n");
  fprintf(out, "\\legend{D,F,B,prepsi,total}\n");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void fastsum_print_output_histo_PreRfNfT(FILE *out, s_testset testset)
{
  int i, size = testset.nresults;
  char hostname[1025];
  char plottitle[1025];

#ifdef HAVE_GETHOSTNAME
  if (gethostname(hostname, 1024) != 0)
#endif
    strncpy(hostname, "unnamed", 1024);

  fastsum_get_plot_title_minus_indep(plottitle, 1024, hostname, testset.param,
      0);

  fprintf(out, "\\begin{tikzpicture}\n");
  fprintf(out, "\\begin{axis}[");
  fprintf(out, "width=0.9\\textwidth, height=0.6\\textwidth, ");
  fprintf(out, "symbolic x coords={");
  for (i = 0; i < size; i++)
    if (i > 0)
      fprintf(out, ",%d", testset.results[i].nthreads);
    else
      fprintf(out, "%d", testset.results[i].nthreads);

  fprintf(out,
      "}, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Time in s, xtick=data, legend style={legend columns=1}, ybar, bar width=7pt, ymajorgrids=true, yminorgrids=true, minor y tick num=1, ");
  fprintf(out, " title={%s}", plottitle);
  fprintf(out, " ]\n");
  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[1].avg + testset.results[i].resval[2].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[3].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[4].avg + testset.results[i].resval[5].avg
            + testset.results[i].resval[6].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[7].avg);
  fprintf(out, "};\n");

  fprintf(out, "\\addplot coordinates {");
  for (i = 0; i < size; i++)
    fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
        testset.results[i].resval[9].avg - testset.results[i].resval[0].avg);
  fprintf(out, "};\n");
  fprintf(out,
      "\\legend{prepsi (step 1b),init nearfield (step 1c),far field (steps 2a-c),nearfield (step 2d),total $-$ step 1a}\n");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void fastsum_print_output_speedup_total_minus_indep(FILE *out,
    s_testset *testsets, int ntestsets)
{
  int i, t;
  char hostname[1025];
  char plottitle[1025];
  unsigned int diff_mask = fastsum_determine_different_parameters(testsets,
      ntestsets);

#ifdef HAVE_GETHOSTNAME
  if (gethostname(hostname, 1024) != 0)
#endif
    strncpy(hostname, "unnamed", 1024);

  fastsum_get_plot_title_minus_indep(plottitle, 1024, hostname,
      testsets[0].param, diff_mask | MASK_FSUM_WINM);

  fprintf(out, "\\begin{tikzpicture}\n");
  fprintf(out, "\\begin{axis}[");
  fprintf(out,
      "width=0.9\\textwidth, height=0.6\\textwidth, x tick label style={ /pgf/number format/1000 sep=}, xlabel=Number of threads, ylabel=Speedup, xtick=data, legend style={ legend pos = north west, legend columns=1}, ymajorgrids=true, yminorgrids=true, minor y tick num=4, ");
  fprintf(out, " title={%s}", plottitle);
  fprintf(out, " ]\n");

  for (t = 0; t < ntestsets; t++)
  {
    s_testset testset = testsets[t];

    R tref = K(0.0);
    for (i = 0; i < testset.nresults; i++)
      if (testset.results[i].nthreads == 1)
        tref = testset.results[i].resval[9].avg
            - testset.results[i].resval[0].avg;

    fprintf(out, "\\addplot coordinates {");
    for (i = 0; i < testset.nresults; i++)
      fprintf(out, "(%d, %.6" __FES__ ") ", testset.results[i].nthreads,
          tref
              / (testset.results[i].resval[9].avg
                  - testset.results[i].resval[0].avg));
    fprintf(out, "};\n");

    for (i = 0; i < testset.nresults; i++)
    {
      fprintf(stderr, "%d:%.3" __FIS__ "  ", testset.results[i].nthreads,
          tref
              / (testset.results[i].resval[9].avg
                  - testset.results[i].resval[0].avg));
    }
    fprintf(stderr, "\n\n");
  }

  fprintf(out, "\\legend{{");
  for (t = 0; t < ntestsets; t++)
  {
    char title[256];
    if (t > 0)
      fprintf(out, "},{");
    fastsum_get_plot_title_minus_indep(title, 255, "", testsets[t].param,
        ~(diff_mask | MASK_FSUM_WINM));
    fprintf(out, "%s", title);
  }
  fprintf(out, "}}\n");
  fprintf(out, "\\end{axis}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fprintf(out, "\n\n");

  fflush(out);
}

void run_testset(s_testset *testset, int d, int L, int M, int n, int m, int p,
    char *kernel_name, R c, R eps_I, R eps_B,
    int *nthreads_array, int n_threads_array_size)
{
  int i;
  testset->param.d = d;
  testset->param.L = L;
  testset->param.M = M;
  testset->param.n = n;
  testset->param.m = m;
  testset->param.p = p;
  testset->param.kernel_name = kernel_name;
  testset->param.c = c;
  testset->param.eps_I = eps_I;
  testset->param.eps_B = eps_B;

  testset->results = (s_result*) NFFT(malloc)(
      (size_t)(n_threads_array_size) * sizeof(s_result));
  testset->nresults = n_threads_array_size;

  run_test_create(testset->param.d, testset->param.L, testset->param.M);
  for (i = 0; i < n_threads_array_size; i++)
  {
    testset->results[i].nthreads = nthreads_array[i];
    run_test(testset->results[i].resval, NREPEAT, testset->param.n,
        testset->param.m, testset->param.p, testset->param.kernel_name,
        testset->param.c, testset->param.eps_I, testset->param.eps_B,
        testset->results[i].nthreads);
  }

}

void test1(int *nthreads_array, int n_threads_array_size)
{
  s_testset testsets[1];

#if defined MEASURE_TIME && defined MEASURE_TIME_FFTW
  run_testset(&testsets[0], 3, 100000, 100000, 128, 4, 7, "one_over_x", K(0.0), K(0.03125), K(0.03125), nthreads_array, n_threads_array_size);

  fastsum_print_output_speedup_total_minus_indep(file_out_tex, testsets, 1);

  fastsum_print_output_histo_PreRfNfT(file_out_tex, testsets[0]);

  nfft_adjoint_print_output_histo_DFBRT(file_out_tex, testsets[0]);

  nfft_trafo_print_output_histo_DFBRT(file_out_tex, testsets[0]);
#endif
}

int main(int argc, char** argv)
{
  int *nthreads_array;
  int n_threads_array_size = get_nthreads_array(&nthreads_array);
  int k;

#if !(defined MEASURE_TIME && defined MEASURE_TIME_FFTW)
  fprintf(stderr, "WARNING: Detailed time measurements are not activated.\n");
  fprintf(stderr, "Please re-run the configure script with options\n");
  fprintf(stderr,
      "--enable-measure-time --enable-measure-time-fftw --enable-openmp\n");
  fprintf(stderr, "and run \"make clean all\"\n\n");
#endif

  for (k = 0; k < n_threads_array_size; k++)
    fprintf(stderr, "%d ", nthreads_array[k]);
  fprintf(stderr, "\n");

  file_out_tex = fopen("fastsum_benchomp_results_plots.tex", "w");

  test1(nthreads_array, n_threads_array_size);

  fclose(file_out_tex);

  return EXIT_SUCCESS;
}

