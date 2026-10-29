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

/* Tiles of the grid of the matrices B and B^T of NFFT, NFCT and NFST (see
 * infft.h): their layout, the node range of each tile in the sorted node
 * order, and the threaded walk over the tiles in colour phases that B^T uses
 * to spread the nodes of different tiles at the same time. */

#include "config.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "nfft3.h"
#include "infft.h"

INT Y(tile_layout)(const INT d, const INT L, const INT *ext, INT *edge,
    INT *ntiles)
{
  INT t, total = 1;
  /* NFFT.jl's default tile edges (src/binning.jl), at least L-1 so that the
   * windows of two tiles of the same colour are disjoint */
  INT e = d == 1 ? 1024 : d == 2 ? 64 : d == 3 ? (L <= 6 ? 16 : L <= 10 ? 12
      : 10) : L;

  e = MAX(e, L - 1);

  for (t = 0; t < d; t++)
  {
    /* d >= 3: the tiles span the last axis. Cubes spread directly into the
     * grid read and write each grid point ((e+L-1)/e)^d times, 5 times for
     * d = 3 and m = 4, from memory, as the tiles of one colour lie apart;
     * columns cut that to ((e+L-1)/e)^(d-1) with long contiguous rows
     * (measured: 3D adjoint with 4 threads twice as fast as with cubes). */
    edge[t] = (d >= 3 && t == d - 1) ? ext[t] : MAX(MIN(e, ext[t]), 1);
    ntiles[t] = (ext[t] + edge[t] - 1) / edge[t];
    total *= ntiles[t];
  }

  return total;
}

void Y(tile_ptr_from_keys)(const INT M, const INT *ar_x, const INT tilevol,
    const INT nbuckets, INT *tile_ptr)
{
  INT k;

  if (M == 0)
  {
    for (k = 0; k <= nbuckets; k++)
      tile_ptr[k] = 0;
    return;
  }

  /* tile_ptr[b] is the first k whose bucket is b or later: written by the
   * node k that starts a new bucket, for the buckets up to its own */
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    const INT b = ar_x[2 * k] / tilevol;
    INT a = k == 0 ? -1 : ar_x[2 * (k - 1)] / tilevol, i;

    for (i = a + 1; i <= b; i++)
      tile_ptr[i] = k;
  }

  for (k = ar_x[2 * (M - 1)] / tilevol + 1; k <= nbuckets; k++)
    tile_ptr[k] = M;
}

/* Number i of the tiles of colour c (bit t of c: parity of the tile
 * coordinate on axis t), in the order of the tile index; cnt[2t+p]: number of
 * tile coordinates of parity p on axis t. */
static INT tile_of_colour(const INT d, const INT *ntiles, const INT *cnt,
    const INT c, INT i)
{
  INT t, tile = 0, stride = 1;

  for (t = d - 1; t >= 0; t--)
  {
    const INT p = (c >> t) & 1, q = i % cnt[2 * t + p];
    i /= cnt[2 * t + p];
    tile += (2 * q + p) * stride;
    stride *= ntiles[t];
  }

  return tile;
}

void Y(tiles_by_colour)(const INT d, const INT *ntiles, const INT *tile_ptr,
    Y(tile_fn) fn, void *ctx)
{
  const INT ncol = (INT)1 << d;
  /* on the heap: no run-time sized arrays in the OpenMP region (GCC 13 and
   * newer abort on them on aarch64-apple-darwin) */
  INT *cnt = (INT*) Y(malloc)((size_t)(2 * d) * sizeof(INT));
  INT c, t;

  for (t = 0; t < d; t++)
  {
    cnt[2 * t] = (ntiles[t] + 1) / 2;
    cnt[2 * t + 1] = ntiles[t] / 2;
  }

#ifdef _OPENMP
  #pragma omp parallel default(shared) private(c, t)
#endif
  for (c = 0; c < ncol; c++)
  {
    INT n = 1, i;

    for (t = 0; t < d; t++)
      n *= cnt[2 * t + ((c >> t) & 1)];

    /* the implicit barrier at the end separates the colours */
#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 1)
#endif
    for (i = 0; i < n; i++)
    {
      const INT b = tile_of_colour(d, ntiles, cnt, c, i);
      if (tile_ptr[b] < tile_ptr[b + 1])
        fn(ctx, tile_ptr[b], tile_ptr[b + 1]);
    }
  }

  Y(free)(cnt);
}
