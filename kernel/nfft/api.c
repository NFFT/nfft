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

#include <stdio.h>  /* FILE, fopen, fclose, fread, fwrite, fseek, ftell, rewind */
#include <stdlib.h> /* accessed via Y(malloc)/Y(free) */

/* In the NFFT module, X() names the legacy public API. */
#undef X
#define X(name) NFFT(name)

/* Print the plan tree to the given FILE stream. Uses the file-backend
 * printer (Y(printer_create_file)) which flushes its internal buffer on
 * Y(printer_destroy). */
void Y(fprint_plan)(Y(plan_ng) * p, FILE *out) {
  printer *pr;
  Y(nfft_ensure_registered)
  ();
  pr = Y(printer_create_file)(out);
  Y(plan_ng_print)
  (p, pr);
  Y(printer_destroy)
  (pr);
}

/* Export wisdom to string: first count bytes with the counter printer, then
 * Y(malloc) count+1 bytes, fill with the string printer, NUL-terminate,
 * and return the caller-owned string. */
char *Y(export_wisdom_to_string)(void) {
  planner *pl;
  printer *pr;
  size_t cnt;
  char *s;

  Y(nfft_ensure_registered)
  ();
  pl = Y(the_planner)();

  /* Pass 1: count characters. */
  pr = Y(printer_create_cnt)(&cnt);
  Y(planner_export)
  (pl, pr);
  Y(printer_destroy)
  (pr);

  /* Allocate: cnt characters + NUL. */
  s = (char *)Y(malloc)(cnt + 1);
  if (!s)
    return 0;

  /* Pass 2: fill the buffer. */
  pr = Y(printer_create_str)(s);
  Y(planner_export)
  (pl, pr);
  Y(printer_destroy)
  (pr);
  /* printer_create_str NUL-terminates incrementally; enforce at end. */
  s[cnt] = '\0';

  return s;
}

/* Import wisdom from a NUL-terminated string. Returns 1 on success,
 * 0 on failure (malformed input, signature mismatch, unknown solver).
 * Y(planner_import) restores its snapshot on failure. */
int Y(import_wisdom_from_string)(const char *s) {
  planner *pl;
  scanner *sc;
  int ret;

  Y(nfft_ensure_registered)
  ();
  pl = Y(the_planner)();

  sc = Y(scanner_create_str)(s);
  ret = Y(planner_import)(pl, sc);
  Y(scanner_destroy)
  (sc);

  return ret;
}

/* Export wisdom to a named file. Returns 1 on success, 0 on any fopen
 * or fwrite failure. The FILE* is always closed before returning.
 * Two-pass (count then write via string buffer) so no file printer
 * failure handling is needed beyond fwrite. */
int Y(export_wisdom_to_filename)(const char *filename) {
  planner *pl;
  printer *pr;
  FILE *f;
  size_t cnt;
  char *s;
  int ok;

  Y(nfft_ensure_registered)
  ();
  pl = Y(the_planner)();

  f = fopen(filename, "wb");
  if (!f)
    return 0;

  /* Pass 1: count characters. */
  pr = Y(printer_create_cnt)(&cnt);
  Y(planner_export)
  (pl, pr);
  Y(printer_destroy)
  (pr);

  s = (char *)Y(malloc)(cnt + 1);
  if (!s) {
    fclose(f);
    return 0;
  }

  /* Pass 2: fill the string buffer. */
  pr = Y(printer_create_str)(s);
  Y(planner_export)
  (pl, pr);
  Y(printer_destroy)
  (pr);
  s[cnt] = '\0';

  ok = (fwrite(s, 1, cnt, f) == cnt) ? 1 : 0;
  Y(free)
  (s);
  fclose(f);

  return ok;
}

/* Import wisdom from a named file. Reads the file into a malloc'd buffer
 * and delegates to the string-scanner path (consistent with the string
 * variant's atomicity semantics). Returns 1 on success, 0 on any
 * fopen/fread/parse failure. */
int Y(import_wisdom_from_filename)(const char *filename) {
  FILE *f;
  long len;
  char *buf;
  size_t nread;
  int ret;

  Y(nfft_ensure_registered)
  ();

  f = fopen(filename, "rb");
  if (!f)
    return 0;

  /* Determine file length. */
  if (fseek(f, 0, SEEK_END) != 0) {
    fclose(f);
    return 0;
  }
  len = ftell(f);
  if (len < 0) {
    fclose(f);
    return 0;
  }
  rewind(f);

  buf = (char *)Y(malloc)((size_t)len + 1);
  if (!buf) {
    fclose(f);
    return 0;
  }

  nread = fread(buf, 1, (size_t)len, f);
  fclose(f);

  if ((long)nread != len) {
    Y(free)
    (buf);
    return 0;
  }
  buf[len] = '\0';

  ret = Y(import_wisdom_from_string)(buf);
  Y(free)
  (buf);

  return ret;
}

/* Erase all wisdom (blessed and unblessed) from the process-global
 * planner. */
void Y(forget_wisdom)(void) {
  Y(nfft_ensure_registered)
  ();
  Y(planner_forget)
  (Y(the_planner)(), PLNR_FORGET_ALL);
}

/* Public per-process planning timelimit. Negative = unlimited. */
void Y(set_timelimit)(double seconds) {
  Y(planner_set_timelimit)
  (Y(the_planner)(), seconds);
}
