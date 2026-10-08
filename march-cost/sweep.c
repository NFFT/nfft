/* Relative l2 error of 1-D NFFT direct and fast transforms against a long double
 * reference, over N = 2^4 .. 2^maxlog. Build without -ffast-math (Kahan sums).
 * Usage: sweep <maxlog> <M>. Output: kind dir N err. */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nfft3.h"

#ifdef USE_FLOAT
#define P(name) nfftf_##name
typedef float R;
#define PREC "float"
#else
#define P(name) nfft_##name
typedef double R;
#define PREC "double"
#endif

static const long double K2PI = 6.283185307179586476925286766559005768L;
static unsigned long long seed;

static double rnd(void) /* uniform in [-0.5, 0.5), same sequence for every build */
{
  seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
  return (double)(seed >> 11) / 9007199254740992.0 - 0.5;
}

typedef struct { long double s, c; } kahan;
static void kadd(kahan *k, long double v)
{
  long double y = v - k->c, t = k->s + y;
  k->c = (t - k->s) - y;
  k->s = t;
}

/* e^{sign 2 pi i k x}, phase reduced exactly mod 1 */
static void phase(long double k, long double x, int sign, long double *re, long double *im)
{
  long double r = fmal(k, x, -rintl(k * x));
  *re = cosl(K2PI * r);
  *im = sign * sinl(K2PI * r);
}

/* relative l2 error of n complex values (R[2] layout) against long double ref */
static double relerr(const R (*v)[2], const long double (*ref)[2], int n)
{
  long double num = 0, den = 0;
  for (int i = 0; i < n; i++) {
    long double dr = v[i][0] - ref[i][0], di = v[i][1] - ref[i][1];
    num += dr * dr + di * di;
    den += ref[i][0] * ref[i][0] + ref[i][1] * ref[i][1];
  }
  return (double)sqrtl(num / den);
}

int main(int argc, char **argv)
{
  int maxlog = argc > 1 ? atoi(argv[1]) : 16, M = argc > 2 ? atoi(argv[2]) : 128;
  for (int lg = 4; lg <= maxlog; lg++) {
    int N = 1 << lg;
    P(plan) p;
    P(init_1d)(&p, N, M);
    seed = 12345 + lg;
    for (int j = 0; j < M; j++)
      p.x[j] = (R)rnd();
    if (p.flags & PRE_ONE_PSI)
      P(precompute_one_psi)(&p);
    R (*fhat)[2] = (R(*)[2])p.f_hat, (*f)[2] = (R(*)[2])p.f;
    R (*fhat0)[2] = malloc(sizeof(R[2]) * N), (*f0)[2] = malloc(sizeof(R[2]) * M);
    long double (*ref_f)[2] = malloc(sizeof(long double[2]) * M);
    long double (*ref_fhat)[2] = malloc(sizeof(long double[2]) * N);
    for (int k = 0; k < N; k++)
      fhat0[k][0] = (R)rnd(), fhat0[k][1] = (R)rnd();
    for (int j = 0; j < M; j++)
      f0[j][0] = (R)rnd(), f0[j][1] = (R)rnd();

    /* trafo: f_j = sum_k fhat_k e^{-2 pi i k x_j} */
    for (int j = 0; j < M; j++) {
      kahan sr = {0, 0}, si = {0, 0};
      for (int k = 0; k < N; k++) {
        long double c, s;
        phase(k - N / 2, p.x[j], -1, &c, &s);
        kadd(&sr, fhat0[k][0] * c - fhat0[k][1] * s);
        kadd(&si, fhat0[k][0] * s + fhat0[k][1] * c);
      }
      ref_f[j][0] = sr.s, ref_f[j][1] = si.s;
    }
    /* adjoint: fhat_k = sum_j f_j e^{+2 pi i k x_j} */
    for (int k = 0; k < N; k++) {
      kahan sr = {0, 0}, si = {0, 0};
      for (int j = 0; j < M; j++) {
        long double c, s;
        phase(k - N / 2, p.x[j], 1, &c, &s);
        kadd(&sr, f0[j][0] * c - f0[j][1] * s);
        kadd(&si, f0[j][0] * s + f0[j][1] * c);
      }
      ref_fhat[k][0] = sr.s, ref_fhat[k][1] = si.s;
    }

    for (int fast = 0; fast < 2; fast++) {
      for (int k = 0; k < N; k++)
        fhat[k][0] = fhat0[k][0], fhat[k][1] = fhat0[k][1];
      if (fast) P(trafo)(&p); else P(trafo_direct)(&p);
      printf("%s %s trafo %d %.6e\n", PREC, fast ? "fast" : "direct", N, relerr(f, ref_f, M));
      for (int j = 0; j < M; j++)
        f[j][0] = f0[j][0], f[j][1] = f0[j][1];
      if (fast) P(adjoint)(&p); else P(adjoint_direct)(&p);
      printf("%s %s adjoint %d %.6e\n", PREC, fast ? "fast" : "direct", N, relerr(fhat, ref_fhat, N));
    }
    fflush(stdout);
    free(fhat0), free(f0), free(ref_f), free(ref_fhat);
    P(finalize)(&p);
  }
  return 0;
}
