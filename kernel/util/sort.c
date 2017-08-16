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


#ifdef _OPENMP
#include <omp.h>
#endif

#include "infft.h"

#define z_swap(_a_, _b_, _t_) do { (_t_) = (_a_); (_a_) = (_b_); (_b_) = (_t_); } while (0)

/**
 * Auxiliary function of radix sort for node indices.
 *
 * \author Michael Hofmann
 */
static inline void sort_node_indices_sort_bubble(INT n, INT *keys)
{
  INT i, j, ti;

  for (i = 0; i < n; ++i)
  {
    j = i;
    while (j > 0 && keys[2 * j + 0] < keys[2 * (j - 1) + 0])
    {
      z_swap(keys[2 * j + 0], keys[2 * (j - 1) + 0], ti);
      z_swap(keys[2 * j + 1], keys[2 * (j - 1) + 1], ti);
      --j;
    }
  }
}

/**
 * Auxiliary function of radix sort for node indices.
 *
 * \author Michael Hofmann
 */
static inline void sort_node_indices_radix_count(INT n, INT *keys, INT shift, INT mask, INT *counts)
{
  INT i, k;

  for (i = 0; i < n; ++i)
  {
    k = (keys[2 * i + 0] >> shift) & mask;
    ++counts[k];
  }
}

/**
 * Auxiliary function of radix sort for node indices.
 *
 * \author Michael Hofmann
 */
static inline void sort_node_indices_radix_rearrange(INT n, INT *keys_in, INT *keys_out, INT shift, INT mask, INT *displs)
{
  INT i, k;

  for (i = 0; i < n; ++i)
  {
    k = (keys_in[2 * i + 0] >> shift) & mask;
    keys_out[2 * displs[k] + 0] = keys_in[2 * i + 0];
    keys_out[2 * displs[k] + 1] = keys_in[2 * i + 1];
    ++displs[k];
  }
}

#define rwidth 9
#define radix_n (1 << rwidth)

/**
 * Radix sort for node indices with OpenMP support.
 *
 * \author Michael Hofmann
 */
void Y(sort_node_indices_radix_lsdf)(INT n, INT *keys0, INT *keys1, INT rhigh)
{
  const INT radix_mask = radix_n - 1;
  const INT rhigh_in = rhigh;

  const INT tmax =
#ifdef _OPENMP
    omp_get_max_threads();
#else
    1;
#endif

  INT *from, *to, *tmp;

  INT i, k, l, h;
  INT *lcounts;

  INT tid = 0, tnum = 1;

  STACK_MALLOC(INT*, lcounts, (size_t)(tmax * radix_n) * sizeof(INT));

  from = keys0;
  to = keys1;

  while (rhigh >= 0)
  {
#ifdef _OPENMP
    #pragma omp parallel private(tid, tnum, i, l, h)
    {
      tid = omp_get_thread_num();
      tnum = omp_get_num_threads();
#endif

      for (i = 0; i < radix_n; ++i) lcounts[tid * radix_n + i] = 0;

      l = (tid * n) / tnum;
      h = ((tid + 1) * n) / tnum;

      sort_node_indices_radix_count(h - l, from + (2 * l), rhigh_in - rhigh, radix_mask, &lcounts[tid * radix_n]);
#ifdef _OPENMP
    }
#endif

    k = 0;
    for (i = 0; i < radix_n; ++i)
    {
      for (l = 0; l < tmax; ++l) lcounts[l * radix_n + i] = (k += lcounts[l * radix_n + i]) - lcounts[l * radix_n + i];
    }

#ifdef _OPENMP
    #pragma omp parallel private(tid, tnum, i, l, h)
    {
      tid = omp_get_thread_num();
      tnum = omp_get_num_threads();
#endif

      l = (tid * n) / tnum;
      h = ((tid + 1) * n) / tnum;

      sort_node_indices_radix_rearrange(h - l, from + (2 * l), to, rhigh_in - rhigh, radix_mask, &lcounts[tid * radix_n]);
#ifdef _OPENMP
    }
#endif

/*    print_keys(n, to);*/

    tmp = from;
    from = to;
    to = tmp;

    rhigh -= rwidth;
  }

  if (to == keys0) memcpy(to, from, (size_t)(n) * 2 * sizeof(INT));

  STACK_FREE(lcounts);
}

/**
 * Radix sort for node indices with OpenMP support.
 *
 * \author Michael Hofmann
 */
void Y(sort_node_indices_radix_msdf)(INT n, INT *keys0, INT *keys1, INT rhigh)
{
  const INT radix_mask = radix_n - 1;

  const INT tmax =
#ifdef _OPENMP
    omp_get_max_threads();
#else
    1;
#endif

  INT i, k, l, h;
  INT *lcounts;

  INT counts[radix_n], displs[radix_n];

  INT tid = 0, tnum = 1;

  STACK_MALLOC(INT*, lcounts, (size_t)(tmax * radix_n) * sizeof(INT));

  rhigh -= rwidth;

#ifdef _OPENMP
  #pragma omp parallel private(tid, tnum, i, l, h)
  {
    tid = omp_get_thread_num();
    tnum = omp_get_num_threads();
#endif

    for (i = 0; i < radix_n; ++i) lcounts[tid * radix_n + i] = 0;

    l = (tid * n) / tnum;
    h = ((tid + 1) * n) / tnum;

    sort_node_indices_radix_count(h - l, keys0 + (2 * l), rhigh + 1, radix_mask, &lcounts[tid * radix_n]);
#ifdef _OPENMP
  }
#endif

  k = 0;
  for (i = 0; i < radix_n; ++i)
  {
    for (l = 0; l < tmax; ++l) lcounts[l * radix_n + i] = (k += lcounts[l * radix_n + i]) - lcounts[l * radix_n + i];

    displs[i] = lcounts[0 * radix_n + i];
    if (i > 0) counts[i - 1] = displs[i] - displs[i - 1];
  }
  counts[radix_n - 1] = n - displs[radix_n - 1];

#ifdef _OPENMP
  #pragma omp parallel private(tid, tnum, i, l, h)
  {
    tid = omp_get_thread_num();
    tnum = omp_get_num_threads();
#endif

    l = (tid * n) / tnum;
    h = ((tid + 1) * n) / tnum;

    sort_node_indices_radix_rearrange(h - l, keys0 + (2 * l), keys1, rhigh + 1, radix_mask, &lcounts[tid * radix_n]);
#ifdef _OPENMP
  }
#endif

  memcpy(keys0, keys1, (size_t)(n) * 2 * sizeof(INT));

  if (rhigh >= 0)
  {
    for (i = 0; i < radix_n; ++i)
    {
      if (counts[i] > 1)
      {
        if (counts[i] > 256)
          Y(sort_node_indices_radix_msdf)(counts[i], keys0 + 2 * displs[i], keys1 + 2 * displs[i], rhigh);
        else
          sort_node_indices_sort_bubble(counts[i], keys0 + 2 * displs[i]);
      }
    }
  }

  STACK_FREE(lcounts);
}
