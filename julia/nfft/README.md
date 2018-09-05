Julia NFFT Interface
====================

General Remarks
---------------

* This interface was developed for Julia 1.0.0 and does not work with previous versions due to significant changes.
* You can enable compilation of the required library by adding `--enable.julia` to the configure call. Please consult the main README for further information.
* The C library does not depend on your Julia installation.

Path
-----

The NFFT module contains an absolute path to the library. You don't need to change anything if you always keep the .so file (or a symbolic link to it) in the same directory as the NFFT module NFFT.jl, otherwise you have to change the lib_path variable accordingly. In this case, you need to watch out for the `pwd()` function since it is executed within the directory of the NFFT.jl module and not the directory where you are currently running Julia.

In order to load the module, you can push the module directory to the load path with `push!`(see simple test files). 
