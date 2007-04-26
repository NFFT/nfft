/* $Id$
 *
 * Copyright (c) 2007 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "config.h"

#if defined(__GNUC__) && (defined(__powerpc__) || defined(__ppc__) || defined(_POWER))
/* The obvious expression a * b + c does not work.  If both x = a * b
   + c and y = a * b - c appear in the source, gcc computes t = a * b,
   x = t + c, y = t - c, thus destroying the fma.

   This peculiar coding seems to do the right thing on all of
   gcc-2.95, gcc-3.1, gcc-3.2, and gcc-3.3.  It does the right thing
   on gcc-3.4 -fno-web (because the ``web'' pass splits the variable
   `x' for the single-assignment form).

   However, gcc-4.0 is a formidable adversary which succeeds in
   pessimizing two fma's into one multiplication and two additions.
   It does it very early in the game---before the optimization passes
   even start.  The only real workaround seems to use fake inline asm
   such as

     asm ("# confuse gcc %0" : "=f"(a) : "0"(a));
     return a * b + c;

   in each of the FMA, FMS, FNMA, and FNMS functions.  However, this
   does not solve the problem either, because two equal asm statements
   count as a common subexpression!  One must use *different* fake asm
   statements:

   in FMA:
     asm ("# confuse gcc for fma %0" : "=f"(a) : "0"(a));

   in FMS:
     asm ("# confuse gcc for fms %0" : "=f"(a) : "0"(a));

   etc.

   After these changes, gcc recalcitrantly generates the fma that was
   in the source to begin with.  However, the extra asm() cruft
   confuses other passes of gcc, notably the instruction scheduler.
   (Of course, one could also generate the fma directly via inline
   asm, but this confuses the scheduler even more.)

   Steven and I have submitted more than one bug report to the gcc
   mailing list over the past few years, to no effect.  Thus, I give
   up.  gcc-4.0 can go to hell.  I'll wait at least until gcc-4.3 is
   out before touching this crap again.
*/
#define E double
static __inline__ E FMA(E a, E b, E c)
{
     E x = a * b;
     x = x + c;
     return x;
}

static __inline__ E FMS(E a, E b, E c)
{
     E x = a * b;
     x = x - c;
     return x;
}

static __inline__ E FNMA(E a, E b, E c)
{
     E x = a * b;
     x = - (x + c);
     return x;
}

static __inline__ E FNMS(E a, E b, E c)
{
     E x = a * b;
     x = - (x - c);
     return x;
}
#else
#define FMA(a, b, c) (((a) * (b)) + (c))
#define FMS(a, b, c) (((a) * (b)) - (c))
#define FNMA(a, b, c) (- (((a) * (b)) + (c)))
#define FNMS(a, b, c) ((c) - ((a) * (b)))
#endif
