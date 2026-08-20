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

/* File backend: buffered, flushed at cleanup. */
#define FILE_PRINTER_BUFSIZE 256

typedef struct
{
  printer base; /* must be first */
  FILE *f;
  char buf[FILE_PRINTER_BUFSIZE];
  int pos;
} printer_file_s;

static void putchr_file(printer *p, char c) {
  printer_file_s *pf = (printer_file_s *)p;
  pf->buf[pf->pos++] = c;
  if (pf->pos == FILE_PRINTER_BUFSIZE) {
    fwrite(pf->buf, 1, (size_t)FILE_PRINTER_BUFSIZE, pf->f);
    pf->pos = 0;
  }
}

static void cleanup_file(printer *p) {
  printer_file_s *pf = (printer_file_s *)p;
  if (pf->pos > 0) {
    fwrite(pf->buf, 1, (size_t)pf->pos, pf->f);
    pf->pos = 0;
  }
}

printer *Y(printer_create_file)(FILE *f) {
  printer_file_s *pf = (printer_file_s *)Y(printer_create)(
      sizeof(printer_file_s), putchr_file, cleanup_file);
  pf->f = f;
  pf->pos = 0;
  return (printer *)pf;
}

/* String backend. buf is borrowed and never grown: the caller must size it for
 * the output plus the NUL, by printing the same thing through the counter
 * backend first. */

typedef struct
{
  printer base; /* must be first */
  char *ptr; /* next write position in the caller's buffer */
} printer_str_s;

static void putchr_str(printer *p, char c) {
  printer_str_s *ps = (printer_str_s *)p;
  *ps->ptr++ = c;
  *ps->ptr = '\0'; /* NUL-terminated after every char */
}

printer *Y(printer_create_str)(char *buf) {
  printer_str_s *ps = (printer_str_s *)Y(printer_create)(
      sizeof(printer_str_s), putchr_str, NULL);
  ps->ptr = buf;
  *buf = '\0';
  return (printer *)ps;
}

/* Counter backend: counts characters without emitting them, excluding the
 * NUL. */

typedef struct
{
  printer base; /* must be first */
  size_t *cnt; /* the caller's counter */
} printer_cnt_s;

static void putchr_cnt(printer *p, char c) {
  printer_cnt_s *pc = (printer_cnt_s *)p;
  UNUSED(c);
  (*pc->cnt)++;
}

printer *Y(printer_create_cnt)(size_t *cnt) {
  printer_cnt_s *pc = (printer_cnt_s *)Y(printer_create)(
      sizeof(printer_cnt_s), putchr_cnt, NULL);
  pc->cnt = cnt;
  *cnt = 0;
  return (printer *)pc;
}
