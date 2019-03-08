Julia NFFT Interface
====================

General Remarks
---------------

* This interface was tested with Julia 1.0.0 and Julia 1.1.0. It does not work with previous versions due to significant changes.
* You can enable compilation of the required library by adding `--enable-julia` to the configure call. Please consult the main README for further information.
* The C library does not depend on your Julia installation.

Path
-----

The NFFT module contains the path to the shared library. You don't have to adjust the variable `lib_path` unless you want to keep the `libnfftjulia.so` in a different directory-

In order to load the NFFT module in Julia, you have to add the directory to the `LOAD_PATH` which would look like `push!(LOAD_PATH, "/path/to/nfft/julia/nfft/")`.
