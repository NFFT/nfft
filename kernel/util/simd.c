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

/* Which of the instruction sets compiled into this build the host can run.
 * See include/simd.h for how the compile-time and the run-time halves fit
 * together. */

#include "infft.h"

/* __builtin_cpu_supports() also checks that the operating system saves the
 * wider register state (XCR0), which is what makes it the right question to
 * ask: a CPU with AVX under a kernel that does not enable it must still take
 * the SSE2 path. */
#if defined(__has_builtin)
#if __has_builtin(__builtin_cpu_supports)
#define NFFT_HAVE_CPU_SUPPORTS 1
#endif
#elif defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8))
#define NFFT_HAVE_CPU_SUPPORTS 1
#endif

int Y(simd_isa_available)(int isa)
{
  switch (isa)
  {
    case NFFT_SIMD_SCALAR:
      return 1;
#if NFFT_SIMD_HAVE_SSE2
    case NFFT_SIMD_SSE2:
#if defined(NFFT_HAVE_CPU_SUPPORTS)
      __builtin_cpu_init();
      return !!__builtin_cpu_supports("sse2");
#else
      /* No run-time query: trust only what the build already targets. */
#if defined(__SSE2__) || defined(__x86_64__)
      return 1;
#else
      return 0;
#endif
#endif
#endif
#if NFFT_SIMD_HAVE_AVX
    case NFFT_SIMD_AVX:
#if defined(NFFT_HAVE_CPU_SUPPORTS)
      __builtin_cpu_init();
      return !!__builtin_cpu_supports("avx");
#elif defined(__AVX__)
      return 1;
#else
      return 0;
#endif
#endif
#if NFFT_SIMD_HAVE_AVX2
    case NFFT_SIMD_AVX2:
#if defined(NFFT_HAVE_CPU_SUPPORTS)
      __builtin_cpu_init();
      return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma") ? 1 : 0;
#elif defined(__AVX2__) && defined(__FMA__)
      return 1;
#else
      return 0;
#endif
#endif
#if NFFT_SIMD_HAVE_NEON
    case NFFT_SIMD_NEON:
      /* Baseline on AArch64. */
      return 1;
#endif
    default:
      return 0;
  }
}

const char *Y(simd_isa_name)(int isa)
{
  switch (isa)
  {
    case NFFT_SIMD_SCALAR:
      return "scalar";
    case NFFT_SIMD_SSE2:
      return "sse2";
    case NFFT_SIMD_AVX:
      return "avx";
    case NFFT_SIMD_AVX2:
      return "avx2";
    case NFFT_SIMD_NEON:
      return "neon";
    default:
      return "unknown";
  }
}

/* Widest variant this build carries that the host also supports. */
static int simd_detect(void)
{
  int isa;

  for (isa = NFFT_SIMD_MAX; isa > NFFT_SIMD_SCALAR; isa--)
    if (Y(simd_isa_available)(isa))
      return isa;

  return NFFT_SIMD_SCALAR;
}

static int simd_from_env(void)
{
  const char *s = getenv("NFFT_SIMD");
  int isa;

  if (s == NULL || *s == '\0')
    return -1;

  if (strcmp(s, "0") == 0 || strcmp(s, "off") == 0 || strcmp(s, "none") == 0)
    return NFFT_SIMD_SCALAR;

  for (isa = NFFT_SIMD_SCALAR; isa <= NFFT_SIMD_NEON; isa++)
    if (strcmp(s, Y(simd_isa_name)(isa)) == 0)
      return isa;

  /* Unrecognised: leave the detection alone rather than silently disable it. */
  return -1;
}

/* Resolved on first use; -1 means "not yet asked". Plans read it while they are
 * being initialised, which is where an application is single-threaded. */
static int simd_active = -1;

int Y(simd_isa)(void)
{
  if (simd_active < 0)
  {
    const int wanted = simd_from_env();

    if (wanted < 0)
      simd_active = simd_detect();
    else
      simd_active = Y(simd_isa_available)(wanted) ? wanted : NFFT_SIMD_SCALAR;
  }

  return simd_active;
}

void Y(simd_force_isa)(int isa)
{
  simd_active = Y(simd_isa_available)(isa) ? isa : NFFT_SIMD_SCALAR;
}
