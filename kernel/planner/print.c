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

static const char hex_digits[] = "0123456789abcdef";

static void print_udec(printer *p, unsigned long v) {
  char d[24];
  int n = 0, i;
  if (v == 0ul) {
    p->putchr(p, '0');
    return;
  }
  while (v != 0ul) {
    d[n++] = (char)('0' + (int)(v % 10ul));
    v /= 10ul;
  }
  for (i = n - 1; i >= 0; i--)
    p->putchr(p, d[i]);
}

static void print_dec(printer *p, long v) {
  unsigned long uv;
  if (v < 0L) {
    p->putchr(p, '-');
    uv = (unsigned long)0ul - (unsigned long)v;
  } else
    uv = (unsigned long)v;
  print_udec(p, uv);
}

/* size_t is the unsigned type of ptrdiff_t's width under C99, so it holds the
 * negation of any INT. */
static void print_INT_dec(printer *p, INT v) {
  size_t uv;
  char d[24];
  int n = 0, i;
  if (v < (INT)0) {
    p->putchr(p, '-');
    uv = (size_t)0 - (size_t)v;
  } else
    uv = (size_t)v;
  if (uv == (size_t)0) {
    p->putchr(p, '0');
    return;
  }
  while (uv != (size_t)0) {
    d[n++] = (char)('0' + (int)(uv % (size_t)10));
    uv /= (size_t)10;
  }
  for (i = n - 1; i >= 0; i--)
    p->putchr(p, d[i]);
}

/* Zero-padded on the left to min_digits. min_digits == 0 with v == 0 emits
 * exactly one '0'. */
static void print_uhex(printer *p, unsigned long v, int min_digits) {
  char d[16];
  int n = 0, i;
  if (v == 0ul && min_digits == 0) {
    p->putchr(p, '0');
    return;
  }
  while (v != 0ul) {
    d[n++] = hex_digits[v & 0xful];
    v >>= 4;
  }
  while (n < min_digits)
    d[n++] = '0';
  for (i = n - 1; i >= 0; i--)
    p->putchr(p, d[i]);
}

/* Directives: %c %s %d (int) %D (INT) %u %x %w (md5uint, 8 hex digits) %p
 * (plan) %P (problem); %( opens an indented line and %) closes it. */
static void vprint(printer *p, const char *fmt, va_list ap) {
  char c;
  int i;
  while ((c = *fmt++) != '\0') {
    if (c != '%') {
      p->putchr(p, c);
      continue;
    }
    c = *fmt++;
    switch (c) {
    case 'c':
      p->putchr(p, (char)va_arg(ap, int));
      break;
    case 's': {
      const char *s = va_arg(ap, const char *);
      if (!s)
        s = "(null)";
      while (*s)
        p->putchr(p, *s++);
      break;
    }
    case 'd':
      print_dec(p, (long)va_arg(ap, int));
      break;
    case 'D':
      print_INT_dec(p, va_arg(ap, INT));
      break;
    case 'u':
      print_udec(p, (unsigned long)va_arg(ap, unsigned));
      break;
    case 'x':
      print_uhex(p, (unsigned long)va_arg(ap, unsigned), 0);
      break;
    case 'w': {
      md5uint w = va_arg(ap, md5uint);
      print_uhex(p, (unsigned long)(w & (md5uint)0xffffffffUL), 8);
      break;
    }
    case 'p': {
      plan *pln = va_arg(ap, plan *);
      if (!pln) {
        const char *s = "(null)";
        while (*s)
          p->putchr(p, *s++);
      } else
        pln->adt->print(pln, p);
      break;
    }
    case 'P': {
      problem *prb = va_arg(ap, problem *);
      if (!prb) {
        const char *s = "(null)";
        while (*s)
          p->putchr(p, *s++);
      } else
        prb->adt->print(prb, p);
      break;
    }
    case '(':
      p->putchr(p, '\n');
      p->indent += p->indent_step;
      for (i = 0; i < p->indent; i++)
        p->putchr(p, ' ');
      break;
    case ')':
      p->indent -= p->indent_step;
      break;
    default:
      A(0); /* unknown directive: programming error */
      break;
    }
  }
}

static void print(printer *p, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  p->vprint(p, fmt, ap);
  va_end(ap);
}

printer *Y(printer_create)(size_t size,
                           void (*putchr)(printer *p, char c), void (*cleanup)(printer *p)) {
  printer *p = (printer *)Y(malloc)(size);
  p->print = print;
  p->vprint = vprint;
  p->putchr = putchr;
  p->cleanup = cleanup;
  p->indent = 0;
  p->indent_step = 2;
  return p;
}

void Y(printer_destroy)(printer *p) {
  if (p->cleanup)
    p->cleanup(p);
  Y(free)
  (p);
}
