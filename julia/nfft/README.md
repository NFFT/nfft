Julia NFFT Interface
====================

General Remarks
---------------

* This interface was developed for Julia 1.0.0 and does not work with previous versions due to significant changes.
* You can enable compilation of the required library by adding `--enable.julia` to the configure call. Please consult the main README for further information.
* The C library does not depend on your Julia installation.

Path
-----

The NFFT module contains an absolute path to the shared library. You have to adjust the variable `lib_path` by replacing `/path/to/nfft/` with the path to your main nfft directory.

In order to load the NFFT module in Julia, you have to add the directory to the `LOAD_PATH` which would look like `push!(LOAD_PATH, "/path/to/nfft/julia/nfft/")`.
