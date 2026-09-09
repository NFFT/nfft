/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* Translation of the caller's request, after FFTW's api/mapflags.c.
 *
 * Patience is the absence of restriction: an impatient request carries more
 * PLNR_NO_* bits and searches a narrower space. The public word names the
 * level; this file turns a level into the restrictions it implies, and into the
 * planner flags handed to the child FFTW plans. Nothing here reads planner or
 * problem state. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* Reserved public bit for a future NFFT_EXHAUSTIVE. The branch below is written
 * against it so that exposing the flag is a one-line change; until then it is
 * unreachable and uncovered, which is a known, accepted cost. */
#define NFFT_EXHAUSTIVE_RESERVED (1U << 6)

/* EXHAUSTIVE implies PATIENT; ESTIMATE denies both. */
static void levels(unsigned planning, unsigned *estimate, unsigned *patient,
                   unsigned *exhaustive)
{
  *estimate = (planning & NFFT_ESTIMATE) ? 1u : 0u;
  *exhaustive = (planning & NFFT_EXHAUSTIVE_RESERVED) ? 1u : 0u;
  *patient = ((planning & NFFT_PATIENT) || *exhaustive) ? 1u : 0u;
  if (*estimate)
    *patient = *exhaustive = 0u;
}

unsigned Y(nfft_map_planning_flags)(unsigned planning)
{
  unsigned estimate, patient, exhaustive;
  unsigned F = 0;

  levels(planning, &estimate, &patient, &exhaustive);

  if (estimate)
    F |= PLNR_ESTIMATE | PLNR_ALLOW_PRUNING;

  if (!patient) /* the fftw2-like block of restrictions */
    F |= PLNR_NO_NONTHREADED | PLNR_BELIEVE_PCOST;

  if (!exhaustive)
    F |= PLNR_NO_UGLY | PLNR_NO_SLOW;

  if (planning & NFFT_NO_NONTHREADED) /* beyond-guru, patience-independent */
    F |= PLNR_NO_NONTHREADED;

  if (planning & NFFT_NO_DIRECT) /* gates are orthogonal to patience */
    F |= PLNR_NO_DIRECT;
  if (planning & NFFT_NO_FAST_NATIVE)
    F |= PLNR_NO_FAST_NATIVE;

  return F;
}

/* The child FFTW plans' flags. Zero means derive from the NFFT patience level,
 * which is what problem_nfft.fftw_flags has always documented. FFTW_PATIENT is
 * the one that does real work: it makes FFTW compare its own threaded and
 * serial candidates instead of preferring the threaded one.
 *
 * FFTW_DESTROY_INPUT is added and FFTW_PRESERVE_INPUT stripped where the child
 * plan is built (kernel/nfft/nfft-nd.c), not here, so those bits stay out of
 * the wisdom key. */
unsigned Y(nfft_derive_fftw_flags)(unsigned planning, unsigned fftw_flags)
{
  unsigned estimate, patient, exhaustive;
  unsigned ff;

  if (fftw_flags != 0u)
    ff = fftw_flags;
  else {
    levels(planning, &estimate, &patient, &exhaustive);
    if (estimate)
      ff = (unsigned)FFTW_ESTIMATE;
    else if (exhaustive)
      ff = (unsigned)FFTW_EXHAUSTIVE;
    else if (patient)
      ff = (unsigned)FFTW_PATIENT;
    else
      ff = (unsigned)FFTW_MEASURE; /* zero */
  }

  /* One source of truth for wisdom-only, on both paths. */
  if (planning & NFFT_WISDOM_ONLY)
    ff |= (unsigned)FFTW_WISDOM_ONLY;
  else
    ff &= ~(unsigned)FFTW_WISDOM_ONLY;

  return ff;
}
