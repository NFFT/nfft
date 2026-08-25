# Make infft.h include its own config.h

## The failure

Configure with Autotools in the source tree, then build with CMake, and the
CMake build silently compiles at the Autotools tree's precision. A float tree
produces a library exporting `nfft_*` and no `nfftf_*`, and nothing warns.

## Why

Two facts combine.

`configure.ac:49` asks for

    AC_CONFIG_HEADERS([include/config.h:include/config.h.in])

and that path is relative to the **build** directory. The flow `CLAUDE.md`
documents is `./configure` from the source root, so build directory equals
source directory and the generated header lands at `<source>/include/config.h`,
beside `infft.h`. An out-of-tree Autotools build puts it under its own
directory and nothing collides -- it is the path that clashes, not the build
directory.

`include/infft.h:25` was

    #include "config.h"

A quoted include searches the directory of the *including* file before any
`-I` path in GCC and Clang. `infft.h` lives in `include/`, so
`include/config.h` won over the `config.h` CMake generated into its own build
directory (`CMakeLists.txt:136`) and put on the `-I` line. No include order
could fix that, because the including file's own directory is consulted first
whatever `-I` says.

## The fix

`#include <config.h>`. An angle include consults only the search path, so the
order on the `-I` line decides, and the build's own directory is put first.

- **CMake**: `kernel/CMakeLists.txt` adds `${NFFT_GENERATED_INCLUDE_DIR}` with
  `BEFORE PRIVATE`, ahead of the public `${PROJECT_SOURCE_DIR}/include`. The
  test targets already listed the generated directory first.
- **Autotools**: nothing to change. Automake already emits
  `DEFAULT_INCLUDES = -I. -I$(top_builddir)/include` and places it ahead of
  `AM_CPPFLAGS` in the compile line, so the build's own header is found first
  in-tree and out-of-tree alike.

`kernel/util/config_guard.c` stays as the compile-time check that the ordering
holds. CMake defines `NFFT_CONFIG_GUARD_PRECISION` to the precision it
configured; the file derives the precision the `config.h` it actually read
selects, and `#error`s if they differ. Builds that leave the macro unset
compile it to nothing, so Autotools is unaffected.

## Verified

With a float `include/config.h` planted in the source tree by
`./configure --enable-float`:

| | result |
|---|---|
| Autotools build of `kernel/util` | compiles, finds its own header through `DEFAULT_INCLUDES` |
| CMake **double** build, `#include <config.h>` | 219 `nfft_*` symbols, **0** `nfftf_*` |
| CMake **double** build, `#include "config.h"` | `config_guard.c` `#error`: shadowed |

The last row is the before state, and it is what this change removes.

## Note

Roughly eighteen other files also write `#include "config.h"`, in
`kernel/*/` and `tests/`. They are harmless: none sits in a directory that
holds a `config.h`, so their quoted include already falls through to the `-I`
path. Only headers under `include/` were exposed, and `infft.h` was the only
one.
