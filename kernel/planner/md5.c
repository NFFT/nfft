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

#include "infft.h"
#include "iplanner.h"

/* Mask to 32 bits: required wherever md5uint is wider than 32 bits. */
#define MD5_M32(x) ((md5uint)(x) & (md5uint)0xffffffffUL)

#define MD5_ROT(x, s) \
  (MD5_M32((md5uint)(x) << (s)) | ((md5uint)(x) >> (32 - (s))))

/* RFC 1321 section 3.4. */
#define MD5_F(b, c, d) (((b) & (c)) | (~(b) & (d)))
#define MD5_G(b, c, d) (((b) & (d)) | ((c) & ~(d)))
#define MD5_H(b, c, d) ((b) ^ (c) ^ (d))
#define MD5_I(b, c, d) ((c) ^ ((b) | ~(d)))

/* a = b + ROT((a + AUX(b,c,d) + x + t), s). The inner mask applies before the
 * rotation, the outer to the final addition. */
#define MD5_FF(a, b, c, d, x, s, t) \
  ((a) = MD5_M32((b) + MD5_ROT(MD5_M32((a) + MD5_F((b), (c), (d)) + (x) + (t)), (s))))
#define MD5_GG(a, b, c, d, x, s, t) \
  ((a) = MD5_M32((b) + MD5_ROT(MD5_M32((a) + MD5_G((b), (c), (d)) + (x) + (t)), (s))))
#define MD5_HH(a, b, c, d, x, s, t) \
  ((a) = MD5_M32((b) + MD5_ROT(MD5_M32((a) + MD5_H((b), (c), (d)) + (x) + (t)), (s))))
#define MD5_II(a, b, c, d, x, s, t) \
  ((a) = MD5_M32((b) + MD5_ROT(MD5_M32((a) + MD5_I((b), (c), (d)) + (x) + (t)), (s))))

/* T[i] = floor(2^32 * |sin(i)|), i = 1..64, stored 0-indexed. */
static const md5uint T64[64] = {
    /* Round 1 */
    0xd76aa478UL, 0xe8c7b756UL, 0x242070dbUL, 0xc1bdceeeUL,
    0xf57c0fafUL, 0x4787c62aUL, 0xa8304613UL, 0xfd469501UL,
    0x698098d8UL, 0x8b44f7afUL, 0xffff5bb1UL, 0x895cd7beUL,
    0x6b901122UL, 0xfd987193UL, 0xa679438eUL, 0x49b40821UL,
    /* Round 2 */
    0xf61e2562UL, 0xc040b340UL, 0x265e5a51UL, 0xe9b6c7aaUL,
    0xd62f105dUL, 0x02441453UL, 0xd8a1e681UL, 0xe7d3fbc8UL,
    0x21e1cde6UL, 0xc33707d6UL, 0xf4d50d87UL, 0x455a14edUL,
    0xa9e3e905UL, 0xfcefa3f8UL, 0x676f02d9UL, 0x8d2a4c8aUL,
    /* Round 3 */
    0xfffa3942UL, 0x8771f681UL, 0x6d9d6122UL, 0xfde5380cUL,
    0xa4beea44UL, 0x4bdecfa9UL, 0xf6bb4b60UL, 0xbebfbc70UL,
    0x289b7ec6UL, 0xeaa127faUL, 0xd4ef3085UL, 0x04881d05UL,
    0xd9d4d039UL, 0xe6db99e5UL, 0x1fa27cf8UL, 0xc4ac5665UL,
    /* Round 4 */
    0xf4292244UL, 0x432aff97UL, 0xab9423a7UL, 0xfc93a039UL,
    0x655b59c3UL, 0x8f0ccc92UL, 0xffeff47dUL, 0x85845dd1UL,
    0x6fa87e4fUL, 0xfe2ce6e0UL, 0xa3014314UL, 0x4e0811a1UL,
    0xf7537e82UL, 0xbd3af235UL, 0x2ad7d2bbUL, 0xeb86d391UL};

static void compress(md5 *ctx) {
  md5uint a, b, c, d, aa, bb, cc, dd;
  md5uint x[16];
  int i;

  for (i = 0; i < 16; i++)
    x[i] = (md5uint)ctx->buf[4 * i] | ((md5uint)ctx->buf[4 * i + 1] << 8) | ((md5uint)ctx->buf[4 * i + 2] << 16) | ((md5uint)ctx->buf[4 * i + 3] << 24);

  aa = a = ctx->s[0];
  bb = b = ctx->s[1];
  cc = c = ctx->s[2];
  dd = d = ctx->s[3];

  /* Round 1 */
  MD5_FF(a, b, c, d, x[0], 7, T64[0]);
  MD5_FF(d, a, b, c, x[1], 12, T64[1]);
  MD5_FF(c, d, a, b, x[2], 17, T64[2]);
  MD5_FF(b, c, d, a, x[3], 22, T64[3]);
  MD5_FF(a, b, c, d, x[4], 7, T64[4]);
  MD5_FF(d, a, b, c, x[5], 12, T64[5]);
  MD5_FF(c, d, a, b, x[6], 17, T64[6]);
  MD5_FF(b, c, d, a, x[7], 22, T64[7]);
  MD5_FF(a, b, c, d, x[8], 7, T64[8]);
  MD5_FF(d, a, b, c, x[9], 12, T64[9]);
  MD5_FF(c, d, a, b, x[10], 17, T64[10]);
  MD5_FF(b, c, d, a, x[11], 22, T64[11]);
  MD5_FF(a, b, c, d, x[12], 7, T64[12]);
  MD5_FF(d, a, b, c, x[13], 12, T64[13]);
  MD5_FF(c, d, a, b, x[14], 17, T64[14]);
  MD5_FF(b, c, d, a, x[15], 22, T64[15]);

  /* Round 2 */
  MD5_GG(a, b, c, d, x[1], 5, T64[16]);
  MD5_GG(d, a, b, c, x[6], 9, T64[17]);
  MD5_GG(c, d, a, b, x[11], 14, T64[18]);
  MD5_GG(b, c, d, a, x[0], 20, T64[19]);
  MD5_GG(a, b, c, d, x[5], 5, T64[20]);
  MD5_GG(d, a, b, c, x[10], 9, T64[21]);
  MD5_GG(c, d, a, b, x[15], 14, T64[22]);
  MD5_GG(b, c, d, a, x[4], 20, T64[23]);
  MD5_GG(a, b, c, d, x[9], 5, T64[24]);
  MD5_GG(d, a, b, c, x[14], 9, T64[25]);
  MD5_GG(c, d, a, b, x[3], 14, T64[26]);
  MD5_GG(b, c, d, a, x[8], 20, T64[27]);
  MD5_GG(a, b, c, d, x[13], 5, T64[28]);
  MD5_GG(d, a, b, c, x[2], 9, T64[29]);
  MD5_GG(c, d, a, b, x[7], 14, T64[30]);
  MD5_GG(b, c, d, a, x[12], 20, T64[31]);

  /* Round 3 */
  MD5_HH(a, b, c, d, x[5], 4, T64[32]);
  MD5_HH(d, a, b, c, x[8], 11, T64[33]);
  MD5_HH(c, d, a, b, x[11], 16, T64[34]);
  MD5_HH(b, c, d, a, x[14], 23, T64[35]);
  MD5_HH(a, b, c, d, x[1], 4, T64[36]);
  MD5_HH(d, a, b, c, x[4], 11, T64[37]);
  MD5_HH(c, d, a, b, x[7], 16, T64[38]);
  MD5_HH(b, c, d, a, x[10], 23, T64[39]);
  MD5_HH(a, b, c, d, x[13], 4, T64[40]);
  MD5_HH(d, a, b, c, x[0], 11, T64[41]);
  MD5_HH(c, d, a, b, x[3], 16, T64[42]);
  MD5_HH(b, c, d, a, x[6], 23, T64[43]);
  MD5_HH(a, b, c, d, x[9], 4, T64[44]);
  MD5_HH(d, a, b, c, x[12], 11, T64[45]);
  MD5_HH(c, d, a, b, x[15], 16, T64[46]);
  MD5_HH(b, c, d, a, x[2], 23, T64[47]);

  /* Round 4 */
  MD5_II(a, b, c, d, x[0], 6, T64[48]);
  MD5_II(d, a, b, c, x[7], 10, T64[49]);
  MD5_II(c, d, a, b, x[14], 15, T64[50]);
  MD5_II(b, c, d, a, x[5], 21, T64[51]);
  MD5_II(a, b, c, d, x[12], 6, T64[52]);
  MD5_II(d, a, b, c, x[3], 10, T64[53]);
  MD5_II(c, d, a, b, x[10], 15, T64[54]);
  MD5_II(b, c, d, a, x[1], 21, T64[55]);
  MD5_II(a, b, c, d, x[8], 6, T64[56]);
  MD5_II(d, a, b, c, x[15], 10, T64[57]);
  MD5_II(c, d, a, b, x[6], 15, T64[58]);
  MD5_II(b, c, d, a, x[13], 21, T64[59]);
  MD5_II(a, b, c, d, x[4], 6, T64[60]);
  MD5_II(d, a, b, c, x[11], 10, T64[61]);
  MD5_II(c, d, a, b, x[2], 15, T64[62]);
  MD5_II(b, c, d, a, x[9], 21, T64[63]);

  /* Add back the saved state. */
  ctx->s[0] = MD5_M32(aa + a);
  ctx->s[1] = MD5_M32(bb + b);
  ctx->s[2] = MD5_M32(cc + c);
  ctx->s[3] = MD5_M32(dd + d);
}

void Y(md5_begin)(md5 *ctx) {
  ctx->s[0] = (md5uint)0x67452301UL;
  ctx->s[1] = (md5uint)0xefcdab89UL;
  ctx->s[2] = (md5uint)0x98badcfeUL;
  ctx->s[3] = (md5uint)0x10325476UL;
  ctx->len = 0;
}

void Y(md5_put_byte)(md5 *ctx, unsigned char b) {
  ctx->buf[ctx->len % 64] = b;
  ctx->len++;
  if (ctx->len % 64 == 0)
    compress(ctx);
}

void Y(md5_put_bytes)(md5 *ctx, const void *data, size_t n) {
  const unsigned char *p = (const unsigned char *)data;
  size_t k;
  for (k = 0; k < n; k++)
    Y(md5_put_byte)
  (ctx, p[k]);
}

void Y(md5_put_str)(md5 *ctx, const char *s) {
  /* The terminating '\0' is fed too, so "ab" and "a","b" hash apart. */
  for (;;) {
    unsigned char b = (unsigned char)*s++;
    Y(md5_put_byte)
    (ctx, b);
    if (b == '\0')
      break;
  }
}

void Y(md5_put_int)(md5 *ctx, int i) {
  Y(md5_put_bytes)
  (ctx, &i, sizeof(i));
}

void Y(md5_put_INT)(md5 *ctx, INT i) {
  Y(md5_put_bytes)
  (ctx, &i, sizeof(i));
}

void Y(md5_put_unsigned)(md5 *ctx, unsigned u) {
  Y(md5_put_bytes)
  (ctx, &u, sizeof(u));
}

void Y(md5_end)(md5 *ctx) {
  unsigned save_len = ctx->len;
  unsigned bits_lo;

  /* The bit count goes out as 64-bit LE with a zero high word: only the low 32
   * bits are computed. save_len * 8u overflows past ~512 MB of input, which
   * would silently diverge from RFC 1321. */
  A(ctx->len < 0x20000000u);

  /* Append bit '1'. */
  Y(md5_put_byte)
  (ctx, (unsigned char)0x80);

  /* Pad to 56 (mod 64), leaving 8 bytes for the length. */
  while (ctx->len % 64 != 56)
    Y(md5_put_byte)
  (ctx, (unsigned char)0x00);

  bits_lo = save_len * 8u;
  Y(md5_put_byte)
  (ctx, (unsigned char)(bits_lo & 0xffu));
  Y(md5_put_byte)
  (ctx, (unsigned char)((bits_lo >> 8) & 0xffu));
  Y(md5_put_byte)
  (ctx, (unsigned char)((bits_lo >> 16) & 0xffu));
  Y(md5_put_byte)
  (ctx, (unsigned char)((bits_lo >> 24) & 0xffu));
  Y(md5_put_byte)
  (ctx, (unsigned char)0x00);
  Y(md5_put_byte)
  (ctx, (unsigned char)0x00);
  Y(md5_put_byte)
  (ctx, (unsigned char)0x00);
  Y(md5_put_byte)
  (ctx, (unsigned char)0x00);
}
