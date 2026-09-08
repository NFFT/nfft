# Code Conventions Used Internally by NFFT3 (not in API)

## Headers

Header            | Installed | Contents                                     | Who includes it
:-----------------|:----------|:---------------------------------------------|:---------------
`nfft3.h`         | yes       | transform APIs, `nfft_malloc`/`nfft_free`    | anyone
`nfft3mp.h`       | yes       | precision-agnostic aliases: `NFFT_R`, `NFFT_C`, `NFFT_K(...)`, `NFFT(...)`, `FFTW(...)`, format strings | anyone writing precision-agnostic code
`nfft3util.h`     | yes       | helpers that belong to no transform: random data, printing, timing, error norms, thread count | anyone
`infft.h`         | no        | the internal names below (`R`, `C`, `Y`, `X`, `A`, `CK`, `DM`, window macros, `ticks`) | `kernel/` and `tests/` only

Examples, applications and the Julia and Matlab bindings use the three
installed headers. They must not include `infft.h`.

The rest of this document describes `infft.h`'s names.

### Include order

`<complex.h>` must be visible before anything pulls in `<fftw3.h>`, or
`fftw_complex`, and therefore `NFFT_C`, is a two-element array of reals that no
arithmetic works on. This is FFTW's contract, not ours. Canonical order:

```c
#include <complex.h>
#include <nfft3.h>
#define NFFT_PRECISION_SINGLE   /* or -DNFFT_PRECISION_* on the command line */
#include <nfft3mp.h>
#include <nfft3util.h>
```

Autotools and CMake both pass the precision macro as `-D`, so in-tree sources
skip the `#define` step.

`nfft3mp.h` includes nothing of its own. `NFFT_C` is a macro that expands to an
`fftw_complex` spelling, so `<fftw3.h>` must already be visible; `nfft3.h` is
what pulls it in. `nfft3mp.h` stops with an `#error` if it is not.

A file that uses a `config.h` macro (`HAVE_FFTW_THREADS`, `MEASURE_TIME`,
`GAUSSIAN`, ...) must `#include "config.h"` itself. The public headers do not
pull it in; only `infft.h` does.



## Common Names

Symbol      | Meaning / Role
:-----------|:--------------------------------------------------------
`R`         | real type, (typically `fftw_real`)
`E`         | real type for local variables (possibly extra precision)
`C`         | complex type
`A`         | assert
`CK`        | check
`X(...)`    | used for mangling of external names (see below)
`Y(...)`    | used for mangling of external names (see below)
`FFTW(...)` | used for mangling of FFTW API names



## Name Mangling

Use `Y(foo)` for external names instead of `nfft_foo`.

`Y(foo)` expands to `nfftf_foo`, `nfftl_foo`, or `nfft_foo`, depending on the
precision.

`X(foo)` also expands to the corresponding prefix, but is local to each module.
I. e. `X(foo)` will expand to `nfct_foo` in the NFCT module, but `Y(foo)` will
always be `nfft_foo`.

`FFTW(foo)` expands to `fftw_foo`.

Names which are not exported do not need to be mangled.



## Indentation and Spacing

Indentation is processed using two spaces:

```
// Level 1.
  // Level 2.
```

Operators and operands, as well, should be padded with a single space character
such that every token is separated from neighbouring ones by a space:

```
int a = 0x1 + 0x2;
```


## Overall Code Style

### General Style Convention

In general, the BSD style should be applied for new lines of code:

```
int foo(int x, int y, int z)
{
  if (x < foo(y, z))
  {
    var = bar[4] + 5;
  }
  else
  {
    while (z)
    {
      var += foo(z, z);
      z--;
    }
    return ++x + bar();
  }
}
```

If desired, the style might be combined with the GNU style as follows:

```
int foo (int x, int y, int z)
{
  if (x < foo (y, z))
    var = bar[4] + 5;
  else
  {
    while (z)
    {
      var += foo (z, z);
      z--;
    }
    return ++x + bar ();
  }
}
```

One style should be applied consequently when writing code.  In case that a
source file was enriched with new lines, these lines should apply the chosen
style.  A refactoring of already existing lines is not required.  Should the
change affect already existing sections of the file, they should be refactored
such that they apply the chosen style.  If a change affects lines which already
apply one of the styles above, the change should follow the convention of the
already existing lines.
