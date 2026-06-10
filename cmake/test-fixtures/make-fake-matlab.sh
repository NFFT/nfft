#!/bin/sh
# Generate the stub MATLAB shared libs for the fake-matlab fixture.
# Usage: make-fake-matlab.sh <fixture-dir>   (defaults to this script's dir/fake-matlab)
set -e
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root="${1:-$here/fake-matlab}"
libdir="$root/bin/glnxa64"
mkdir -p "$libdir"
cc -shared -fPIC -o "$libdir/libmx.so"  "$root/stub.c"
cc -shared -fPIC -o "$libdir/libmex.so" "$root/stub.c"
# also the libs FindMatlab probes for some components (empty is fine)
for L in mat eng MatlabDataArray MatlabEngine; do
  cc -shared -fPIC -o "$libdir/lib$L.so" "$root/stub.c"
done
# libmwfftw3 = the REAL system FFTW (same ABI) so the kernel's fftw_* symbols
# resolve. Name it the VERSIONED libmwfftw3.so.3 (NOT .so) so the fixture exercises
# the finder's versioned-name + build-local-symlink path, exactly as real MATLAB
# ships it. (The copied lib's internal SONAME stays libfftw3.so.3 — fine; the stub
# is a link-time wiring test, not a runtime test.)
realfftw=$(ls /usr/lib/*/libfftw3.so.3 /usr/lib/libfftw3.so.3 2>/dev/null | head -1)
[ -n "$realfftw" ] || { echo "system libfftw3.so.3 not found" >&2; exit 1; }
cp "$realfftw" "$libdir/libmwfftw3.so.3"
# mexext helper (FindMatlab may call it)
printf '#!/bin/sh\necho mexa64\n' > "$root/bin/mexext"; chmod +x "$root/bin/mexext"
echo "stub MATLAB tree ready at $root"
