/* Worst-case error of the 1-D Kaiser-Bessel NFFT, one binary per precision
 * (-DSWEEP_PREC=0,1,2). The sup of max_j |f_j - s_j| / ||f_hat||_1 is attained
 * at f_hat = e_k, so one transform per probed k against the references
 * floorgen.py writes; ||e_k||_1 = 1, so the error is the max over the nodes.
 * PRE_PHI_HUT | PRE_PSI.
 *
 * usage: floor <ref-file> [m_max] [sigmas] [nlo,nhi,m]
 * rows:  prec,sigma_req,sigma,N,n,M,m,err,argk */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nfft3.h"

#if SWEEP_PREC == 0
typedef float R;
typedef float _Complex C;
#define NF(name) nfftf_##name
#define PREC_NAME "float"
#define CRE crealf
#define CIM cimagf
#elif SWEEP_PREC == 1
typedef double R;
typedef double _Complex C;
#define NF(name) nfft_##name
#define PREC_NAME "double"
#define CRE creal
#define CIM cimag
#elif SWEEP_PREC == 2
typedef long double R;
typedef long double _Complex C;
#define NF(name) nfftl_##name
#define PREC_NAME "long_double"
#define CRE creall
#define CIM cimagl
#else
#error set SWEEP_PREC
#endif

typedef long double LR;
typedef long double _Complex LC;

static int M;
static LR *xl;
static R *x;

static void die(const char *s)
{
  fprintf(stderr, "floor: %s\n", s);
  exit(1);
}

/* One N block of the reference file: the probed k and the M references each. */
struct block
{
  int N, nk;
  int *k;
  LC *ref; /* nk * M */
};

static int read_blocks(const char *path, struct block **out)
{
  FILE *f = fopen(path, "r");
  struct block *b = NULL;
  int nb = 0, cap = 0, i, j;
  char tag[16];

  if (!f)
    die("cannot open reference file");
  if (fscanf(f, "%15s %d", tag, &M) != 2 || strcmp(tag, "M"))
    die("bad header");
  xl = (LR *)malloc((size_t)M * sizeof(LR));
  x = (R *)NF(malloc)((size_t)M * sizeof(R));
  if (fscanf(f, "%15s", tag) != 1 || strcmp(tag, "X"))
    die("bad node line");
  for (j = 0; j < M; j++)
  {
    char v[160];
    if (fscanf(f, "%159s", v) != 1)
      die("short node line");
    xl[j] = strtold(v, NULL);
    x[j] = (R)xl[j];
  }

  while (fscanf(f, "%15s", tag) == 1)
  {
    struct block *cur;
    if (strcmp(tag, "N"))
      die("expected N block");
    if (nb == cap)
    {
      cap = cap ? 2 * cap : 8;
      b = (struct block *)realloc(b, (size_t)cap * sizeof(*b));
    }
    cur = &b[nb++];
    if (fscanf(f, "%d %d", &cur->N, &cur->nk) != 2)
      die("bad N block");
    cur->k = (int *)malloc((size_t)cur->nk * sizeof(int));
    cur->ref = (LC *)malloc((size_t)cur->nk * (size_t)M * sizeof(LC));
    if (fscanf(f, "%15s", tag) != 1 || strcmp(tag, "K"))
      die("bad K line");
    for (i = 0; i < cur->nk; i++)
      if (fscanf(f, "%d", &cur->k[i]) != 1)
        die("short K line");
    for (i = 0; i < cur->nk; i++)
      for (j = 0; j < M; j++)
      {
        char re[64], im[64];
        if (fscanf(f, "%63s %63s", re, im) != 2)
          die("short reference");
        cur->ref[(size_t)i * M + j] = strtold(re, NULL) + I * strtold(im, NULL);
      }
  }
  fclose(f);
  *out = b;
  return nb;
}

/* Largest m the fast solver accepts: N > m and n > 2m+2. */
static int m_ok_max(int N, int n)
{
  int a = N - 1, b = (n - 3) / 2;
  return a < b ? a : b;
}

/* max over the probed k and the nodes of |f - s|, which is the measure itself
 * because every probe has ||f_hat||_1 = 1. */
static LR measure(const struct block *blk, int n, int m, int *argk)
{
  NF(plan) p;
  int Ni = blk->N, ni = n, i, j;
  LR worst = 0.0L;

  *argk = 0;
  NF(init_guru)(&p, 1, &Ni, M, &ni, m,
                PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
                    FFTW_INIT | FFT_OUT_OF_PLACE,
                FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  for (j = 0; j < M; j++)
    p.x[j] = x[j];
  if (p.flags & PRE_ONE_PSI)
    NF(precompute_one_psi)(&p);

  for (i = 0; i < blk->nk; i++)
  {
    const LC *ref = blk->ref + (size_t)i * M;
    LR wk = 0.0L;
    memset(p.f_hat, 0, (size_t)blk->N * sizeof(C));
    p.f_hat[blk->k[i] + blk->N / 2] = (C)1.0;
    NF(trafo)(&p);
    for (j = 0; j < M; j++)
    {
      LC a = (LR)CRE(p.f[j]) + I * (LR)CIM(p.f[j]);
      LR e = cabsl(a - ref[j]);
      if (!(e <= wk)) /* NaN propagates: the float window overflows */
        wk = e;
    }
    if (!(wk <= worst))
    {
      worst = wk;
      *argk = blk->k[i];
    }
  }
  NF(finalize)(&p);
  return worst;
}

int main(int argc, char **argv)
{
  double sigmas[32] = {1.5, 2.0, 3.0, 4.0};
  int n_sigma = 4;
  struct block *blk;
  int nb, ib, is, m, i, j;
  const int m_max_cap = (argc > 2 && argv[2][0]) ? atoi(argv[2]) : 40;

  if (argc < 2)
    die("usage: floor <ref-file> [m_max] [sigmas] [nlo,nhi,m]");
  nb = read_blocks(argv[1], &blk);

  if (argc > 3 && argv[3][0])
  {
    char *t = strtok(argv[3], ",");
    for (n_sigma = 0; t && n_sigma < 32; t = strtok(0, ","))
      sigmas[n_sigma++] = atof(t);
  }

  printf("prec,sigma_req,sigma,N,n,M,m,err,argk\n");
  fflush(stdout);

  /* Scan mode: every n in [nlo, nhi] at one cut-off, which is how the grid
   * offset's dependence on n shows up. */
  if (argc > 4 && argv[4][0])
  {
    int nlo = 0, nhi = 0, mfix = 0, n;

    if (sscanf(argv[4], "%d,%d,%d", &nlo, &nhi, &mfix) != 3)
      die("scan wants nlo,nhi,m");
    for (ib = 0; ib < nb; ib++)
      for (n = nlo; n <= nhi; n++)
      {
        int argk;
        LR worst;

        if (n <= blk[ib].N || mfix > m_ok_max(blk[ib].N, n))
          continue;
        worst = measure(&blk[ib], n, mfix, &argk);
        printf("%s,%.4f,%.10f,%d,%d,%d,%d,%.10Le,%d\n", PREC_NAME,
               (double)n / (double)blk[ib].N, (double)n / (double)blk[ib].N,
               blk[ib].N, n, M, mfix, worst, argk);
        fflush(stdout);
      }
    return 0;
  }

  for (ib = 0; ib < nb; ib++)
  {
    const int N = blk[ib].N;

    for (is = 0; is < n_sigma; is++)
    {
      int n = (int)(sigmas[is] * N + 0.5);
      int Ni = N, ni;
      int m_max;
      LR best = 1.0e300L;
      int best_m = 0;

      if (n % 2)
        n++;
      if (n <= N)
        n = N + 2;
      ni = n;
      m_max = m_ok_max(N, n);
      if (m_max > m_max_cap)
        m_max = m_max_cap;

      for (m = 1; m <= m_max; m++)
      {
        LR worst;
        int argk;

        worst = measure(&blk[ib], n, m, &argk);

        printf("%s,%.4f,%.10f,%d,%d,%d,%d,%.10Le,%d\n", PREC_NAME, sigmas[is],
               (double)n / (double)N, N, n, M, m, worst, argk);
        fflush(stdout);

        if (worst < best) /* NaN fails this, so the floor ignores it */
        {
          best = worst;
          best_m = m;
        }
        /* The curve falls to a minimum and then climbs; stop once the climb
         * is unambiguous so the long-double leg stays affordable. */
        if (m >= best_m + 4 && worst > 50.0L * best)
          break;
      }
    }
  }
  return 0;
}
