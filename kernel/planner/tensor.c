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

#include <string.h>

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* tensor / mvdim: rectangular operator geometry for the planner. A tensor is
 * a collection of mvdim factors, each an n_out x n_in matrix over strided
 * input/output vectors. Factor order carries no meaning -- the strides encode
 * the layout -- so tensor_compress canonicalises, making permuted or
 * unit-padded tensors compare and hash equal. */

static INT iabs(INT x) {
  return x < (INT)0 ? -x : x;
}

mvdim Y(mvdim_square)(INT n, INT is, INT os) {
  mvdim d;
  d.n_in = n;
  d.is = is;
  d.n_out = n;
  d.os = os;
  return d;
}

tensor *Y(tensor_create)(int rnk) {
  tensor *t;
  A(rnk >= 0);
  t = (tensor *)Y(malloc)(sizeof(tensor));
  t->rnk = rnk;
  if (rnk > 0)
    t->dims = (mvdim *)Y(malloc)((size_t)rnk * sizeof(mvdim));
  else
    t->dims = 0;
  return t;
}

void Y(tensor_destroy)(tensor *t) {
  if (t) {
    Y(free)
    (t->dims);
    Y(free)
    (t);
  }
}

tensor *Y(tensor_copy)(const tensor *t) {
  tensor *c = Y(tensor_create)(t->rnk);
  if (t->rnk > 0)
    memcpy(c->dims, t->dims, (size_t)t->rnk * sizeof(mvdim));
  return c;
}

int Y(tensor_equal)(const tensor *a, const tensor *b) {
  int i;
  if (a->rnk != b->rnk)
    return 0;
  for (i = 0; i < a->rnk; ++i) {
    const mvdim *x = a->dims + i, *y = b->dims + i;
    if (x->n_in != y->n_in || x->is != y->is || x->n_out != y->n_out || x->os != y->os)
      return 0;
  }
  return 1;
}

int Y(tensor_squarep)(const tensor *t) {
  int i;
  for (i = 0; i < t->rnk; ++i)
    if (t->dims[i].n_in != t->dims[i].n_out)
      return 0;
  return 1;
}

int Y(tensor_kosherp)(const tensor *t) {
  int i;
  if (t->rnk < 0)
    return 0;
  for (i = 0; i < t->rnk; ++i)
    if (t->dims[i].n_in < (INT)1 || t->dims[i].n_out < (INT)1)
      return 0;
  return 1;
}

INT Y(tensor_sz_in)(const tensor *t) {
  INT p = (INT)1;
  int i;
  for (i = 0; i < t->rnk; ++i)
    p *= t->dims[i].n_in;
  return p;
}

INT Y(tensor_sz_out)(const tensor *t) {
  INT p = (INT)1;
  int i;
  for (i = 0; i < t->rnk; ++i)
    p *= t->dims[i].n_out;
  return p;
}

tensor *Y(tensor_adjoint)(const tensor *t) {
  tensor *a = Y(tensor_create)(t->rnk);
  int i;
  for (i = 0; i < t->rnk; ++i) {
    a->dims[i].n_in = t->dims[i].n_out;
    a->dims[i].is = t->dims[i].os;
    a->dims[i].n_out = t->dims[i].n_in;
    a->dims[i].os = t->dims[i].is;
  }
  return a;
}

/* Canonical total order. The signed is/os tiebreakers make it total over
 * fieldwise-distinct factors, so the canonical form does not depend on sort
 * stability. Returns <0 if a precedes b. */
static int dim_cmp(const mvdim *a, const mvdim *b) {
  INT amin, bmin;

  amin = iabs(a->is) < iabs(a->os) ? iabs(a->is) : iabs(a->os);
  bmin = iabs(b->is) < iabs(b->os) ? iabs(b->is) : iabs(b->os);
  if (amin != bmin)
    return amin > bmin ? -1 : 1;

  if (iabs(a->is) != iabs(b->is))
    return iabs(a->is) > iabs(b->is) ? -1 : 1;

  if (iabs(a->os) != iabs(b->os))
    return iabs(a->os) > iabs(b->os) ? -1 : 1;

  if (a->n_in != b->n_in)
    return a->n_in < b->n_in ? -1 : 1;

  if (a->n_out != b->n_out)
    return a->n_out < b->n_out ? -1 : 1;

  if (a->is != b->is)
    return a->is < b->is ? -1 : 1;

  if (a->os != b->os)
    return a->os < b->os ? -1 : 1;

  return 0;
}

static int dim_cmp_qsort(const void *pa, const void *pb) {
  return dim_cmp((const mvdim *)pa, (const mvdim *)pb);
}

static void sort_canonical(mvdim *dims, int n) {
  if (n > 1)
    qsort(dims, (size_t)n, sizeof(mvdim), dim_cmp_qsort);
}

static int unitp(const mvdim *d) {
  return d->n_in == (INT)1 && d->n_out == (INT)1;
}

tensor *Y(tensor_compress)(const tensor *t) {
  tensor *r;
  int i, k = 0;

  for (i = 0; i < t->rnk; ++i)
    if (!unitp(t->dims + i))
      ++k;

  r = Y(tensor_create)(k);
  k = 0;
  for (i = 0; i < t->rnk; ++i)
    if (!unitp(t->dims + i))
      r->dims[k++] = t->dims[i];

  sort_canonical(r->dims, r->rnk);
  return r;
}

static int dim_cmp_is_desc(const void *pa, const void *pb) {
  INT a = iabs(((const mvdim *)pa)->is), b = iabs(((const mvdim *)pb)->is);
  if (a != b)
    return a > b ? -1 : 1;
  return 0;
}

tensor *Y(tensor_compress_contiguous)(const tensor *t) {
  tensor *r;
  mvdim *tmp;
  int i, k = 0;

  A(Y(tensor_kosherp)(t));

  tmp = (mvdim *)Y(malloc)((size_t)(t->rnk > 0 ? t->rnk : 1) * sizeof(mvdim));
  for (i = 0; i < t->rnk; ++i)
    if (!unitp(t->dims + i))
      tmp[k++] = t->dims[i];

  /* Descending |is| puts the outer factor before the inner. When the inner
   * factor has n_in == 1 the two |is| tie and a legal merge can be skipped;
   * the signed nesting test below is exact, so no false merge is possible. */
  if (k > 1)
    qsort(tmp, (size_t)k, sizeof(mvdim), dim_cmp_is_desc);

  /* Merge each adjacent outer/inner pair that nests on both sides. */
  {
    int m = 0; /* write index of the current merged factor */
    for (i = 0; i < k; ++i) {
      if (m > 0) {
        mvdim *a = tmp + (m - 1); /* current accumulator (outer) */
        mvdim *b = tmp + i; /* candidate inner */
        if (a->is == b->is * b->n_in && a->os == b->os * b->n_out) {
          a->n_in = a->n_in * b->n_in;
          a->n_out = a->n_out * b->n_out;
          a->is = b->is;
          a->os = b->os;
          continue;
        }
      }
      tmp[m++] = tmp[i];
    }
    k = m;
  }

  r = Y(tensor_create)(k);
  for (i = 0; i < k; ++i)
    r->dims[i] = tmp[i];
  Y(free)
  (tmp);

  sort_canonical(r->dims, r->rnk);
  return r;
}

/* Feed order: rnk, then per factor in stored order n_in, is, n_out, os. No
 * internal compression -- hash a compressed tensor to get a canonical key. */
void Y(tensor_md5)(md5 *ctx, const tensor *t) {
  int i;
  Y(md5_put_int)
  (ctx, t->rnk);
  for (i = 0; i < t->rnk; ++i) {
    Y(md5_put_INT)
    (ctx, t->dims[i].n_in);
    Y(md5_put_INT)
    (ctx, t->dims[i].is);
    Y(md5_put_INT)
    (ctx, t->dims[i].n_out);
    Y(md5_put_INT)
    (ctx, t->dims[i].os);
  }
}

/* (tensor <rnk> (<n_in> <is> <n_out> <os>) ...); rank 0 prints (tensor 0). */
void Y(tensor_print)(const tensor *t, printer *p) {
  int i;
  p->print(p, "(tensor %d", t->rnk);
  for (i = 0; i < t->rnk; ++i)
    p->print(p, " (%D %D %D %D)", t->dims[i].n_in, t->dims[i].is,
             t->dims[i].n_out, t->dims[i].os);
  p->putchr(p, ')');
}
