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

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* md5 checksum over a node array.  dM is the element count of that very array
 * (rnk*M for a compressed x), never a recomputed d*M. */

void Y(nfft_x_md5)(const R *x, INT dM, md5sig out) {
  md5 ctx;
  Y(md5_begin)
  (&ctx);
  Y(md5_put_bytes)
  (&ctx, x, (size_t)dM * sizeof(R));
  Y(md5_end)
  (&ctx);
  out[0] = ctx.s[0];
  out[1] = ctx.s[1];
  out[2] = ctx.s[2];
  out[3] = ctx.s[3];
}

int Y(nfft_x_verify)(const R *x, INT dM, const md5sig ref) {
  md5sig now;
  Y(nfft_x_md5)
  (x, dM, now);
  return now[0] == ref[0] && now[1] == ref[1] && now[2] == ref[2] &&
         now[3] == ref[3];
}
