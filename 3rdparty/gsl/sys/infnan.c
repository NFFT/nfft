/* sys/infnan.c
 * 
 * Copyright (C) 2001, 2004 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "config.h"
#include <math.h>

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_nan (void)
{
  return gsl_fdiv (0.0, 0.0);
}

double gsl_posinf (void)
{
  return gsl_fdiv (+1.0, 0.0);
}

double gsl_neginf (void)
{
  return gsl_fdiv (-1.0, 0.0);
}


int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

int
gsl_isnan (const double x)
{
  int status = (x != x);
  return status;
}

int
gsl_isinf (const double x)
{
  double y = x - x;
  int s = (y != y);

  if (s && x > 0)
    return +1;
  else if (s && x < 0)
    return -1;
  else
    return 0;
}

int
gsl_finite (const double x)
{
  const double y = x - x;
  int status = (y == y);
  return status;
}
