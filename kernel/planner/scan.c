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

#include <limits.h>

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

static int read_char(scanner *sc) {
  int c;
  if (sc->pushback != EOF) {
    c = sc->pushback;
    sc->pushback = EOF;
    return c;
  }
  return sc->getchr(sc);
}

static void unread_char(scanner *sc, int c) {
  sc->pushback = c;
}

/* Whitespace is any byte <= ' ', so the scan is locale-independent. */
static void skip_ws(scanner *sc) {
  int c;
  for (;;) {
    c = read_char(sc);
    if (c == EOF)
      return;
    if ((unsigned char)c > (unsigned char)' ') {
      unread_char(sc, c);
      return;
    }
  }
}

/* One or more hex digits, either case, into *valp. Returns 1 if at least one
 * digit was consumed. The first non-hex character is pushed back. */
static int read_hex(scanner *sc, unsigned long *valp) {
  unsigned long val = 0;
  int c, n = 0;
  for (;;) {
    c = read_char(sc);
    if (c >= '0' && c <= '9') {
      val = val * 16UL + (unsigned long)(c - '0');
      n++;
    } else if (c >= 'a' && c <= 'f') {
      val = val * 16UL + (unsigned long)(c - 'a' + 10);
      n++;
    } else if (c >= 'A' && c <= 'F') {
      val = val * 16UL + (unsigned long)(c - 'A' + 10);
      n++;
    } else {
      if (c != EOF)
        unread_char(sc, c);
      break;
    }
  }
  *valp = val;
  return (n > 0) ? 1 : 0;
}

/* Directives: %*s (width then string), %d, %x, %w (md5uint). A literal paren
 * skips leading whitespace, a whitespace literal skips a run of any length,
 * every other literal must match exactly. Returns 1 on complete match, 0 on
 * the first mismatch. */
static int vscan(scanner *sc, const char *fmt, va_list ap) {
  char fc;
  while ((fc = *fmt++) != '\0') {
    if (fc == '%') {
      char dc = *fmt++;
      int maxlen = -1;

      if (dc == '*') {
        maxlen = va_arg(ap, int);
        dc = *fmt++;
        A(dc == 's'); /* %* is only defined before %s */
      }

      switch (dc) {
      case 's': {
        /* Up to maxlen non-whitespace, non-paren chars. Leading whitespace
         * is skipped so "(%*s..." works at any input indentation. */
        char *buf = va_arg(ap, char *);
        int c, n = 0;
        A(maxlen >= 0);
        skip_ws(sc);
        for (;;) {
          c = read_char(sc);
          if (c == EOF)
            break;
          if ((unsigned char)c <= (unsigned char)' ' || c == '(' || c == ')') {
            unread_char(sc, c);
            break;
          }
          if (n < maxlen)
            buf[n++] = (char)c;
        }
        buf[n] = '\0';
        break;
      }

      case 'd': {
        /* Optional sign, then one or more decimal digits. */
        int *ip = va_arg(ap, int *);
        int sign = 1, n = 0;
        long long val = 0;
        int c = read_char(sc);
        if (c == '+') {
          sign = 1;
          c = read_char(sc);
        } else if (c == '-') {
          sign = -1;
          c = read_char(sc);
        }
        while (c >= '0' && c <= '9') {
          val = val * 10LL + (long long)(c - '0');
          if (val > (long long)INT_MAX + 1) /* out of range for either sign */
            return 0;
          n++;
          c = read_char(sc);
        }
        if (c != EOF)
          unread_char(sc, c);
        if (n == 0)
          return 0;
        if (sign >= 0 && val > (long long)INT_MAX)
          return 0;
        *ip = (sign >= 0) ? (int)val : (int)(-val); /* negate in long long: no UB */
        break;
      }

      case 'x': {
        unsigned *up = va_arg(ap, unsigned *);
        unsigned long val;
        if (!read_hex(sc, &val))
          return 0;
        *up = (unsigned)val;
        break;
      }

      case 'w': {
        md5uint *mp = va_arg(ap, md5uint *);
        unsigned long val;
        if (!read_hex(sc, &val))
          return 0;
        *mp = (md5uint)(val & (unsigned long)0xffffffffUL);
        break;
      }

      default:
        A(0); /* unknown directive: programming error */
        return 0;
      }
    } else if (fc == '(' || fc == ')') {
      int c;
      skip_ws(sc);
      c = read_char(sc);
      if (c != (int)(unsigned char)fc) {
        if (c != EOF)
          unread_char(sc, c);
        return 0;
      }
    } else if ((unsigned char)fc <= (unsigned char)' ') {
      skip_ws(sc);
    } else {
      int c = read_char(sc);
      if (c != (int)(unsigned char)fc) {
        if (c != EOF)
          unread_char(sc, c);
        return 0;
      }
    }
  }
  return 1;
}

static int scan(scanner *sc, const char *fmt, ...) {
  va_list ap;
  int r;
  va_start(ap, fmt);
  r = sc->vscan(sc, fmt, ap);
  va_end(ap);
  return r;
}

scanner *Y(scanner_create)(size_t size, int (*getchr)(scanner *sc)) {
  scanner *sc = (scanner *)Y(malloc)(size);
  sc->scan = scan;
  sc->vscan = vscan;
  sc->getchr = getchr;
  sc->pushback = EOF;
  return sc;
}

void Y(scanner_destroy)(scanner *sc) {
  Y(free)
  (sc);
}
