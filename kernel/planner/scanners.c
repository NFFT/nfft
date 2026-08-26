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

typedef struct
{
  scanner base; /* must be first */
  FILE *f;
} scanner_file_s;

static int getchr_file(scanner *sc) {
  scanner_file_s *sf = (scanner_file_s *)sc;
  return fgetc(sf->f);
}

scanner *Y(scanner_create_file)(FILE *f) {
  scanner_file_s *sf = (scanner_file_s *)Y(scanner_create)(
      sizeof(scanner_file_s), getchr_file);
  sf->f = f;
  return (scanner *)sf;
}

/* String backend. s is borrowed and must outlive the scanner. */
typedef struct
{
  scanner base; /* must be first */
  const char *ptr; /* current read position */
} scanner_str_s;

static int getchr_str(scanner *sc) {
  scanner_str_s *ss = (scanner_str_s *)sc;
  unsigned char c = (unsigned char)*ss->ptr;
  if (c == '\0')
    return EOF;
  ss->ptr++;
  return (int)c;
}

scanner *Y(scanner_create_str)(const char *s) {
  scanner_str_s *ss = (scanner_str_s *)Y(scanner_create)(
      sizeof(scanner_str_s), getchr_str);
  ss->ptr = s;
  return (scanner *)ss;
}
