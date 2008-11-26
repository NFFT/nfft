#
#  NFFT_LIB_FFTW3
#
#  Check configure options and assign variables related to the fftw3 library.
#
#  If we find the library, set the shell variable `nfft_lib_fftw3' to `yes'.
#  Otherwise, set `nfft_lib_fftw3' to `no'.

AC_DEFUN([NFFT_LIB_FFTW3],
[
  AC_ARG_WITH(fftw3, [AC_HELP_STRING([--with-fftw3=DIR],
  [compile with fftw3 in DIR])], with_fftw3=$withval, with_fftw3="yes")

  AC_ARG_WITH(fftw3-libdir, [AC_HELP_STRING([--with-fftw3-libdir=DIR],
  [compile with fftw3 library directory DIR])], with_fftw3_libdir=$withval, with_fftw3_libdir="yes")

  AC_ARG_WITH(fftw3-includedir, [AC_HELP_STRING([--with-fftw3-includedir=DIR],
  [compile with fftw3 include directory DIR])], with_fftw3_includedir=$withval, with_fftw3_includedir="yes")

  # Save compiler and linker flags.
  saved_LDFLAGS="$LDFLAGS"
  saved_CPPFLAGS="$CPPFLAGS"

  # Check if user specified a valid root directory.
  if test ! -d ${with_fftw3}; then
    fftw3_prefix="/usr/local"
    if test ! -d ${with_fftw3_includedir}; then
      fftw3_includedir="$fftw3_prefix/include"
    else
      fftw3_includedir="$with_fftw3_includedir"
    fi
    if test ! -d ${with_fftw3_libdir}; then
      fftw3_libdir="$fftw3_prefix/lib"
    else
      fftw3_libdir="$with_fftw3_libdir"
    fi
  else
    fftw3_prefix=$with_fftw3
    fftw3_includedir="$fftw3_prefix/include"
    fftw3_libdir="$fftw3_prefix/lib"
  fi

  fftw3_includes="-I$fftw3_includedir"
  fftw3_libs="-L$fftw3_libdir"

  # Modify compiler and linker flags.
  CPPFLAGS="$CPPFLAGS $fftw3_includes"
  LDFLAGS="$LDFLAGS $fftw3_libs"

  # Check if library is present and usable.
  AC_CHECK_HEADER([fftw3.h],
    [AC_CHECK_LIB(fftw3, fftw_execute,
      nfft_lib_fftw3=yes,
      nfft_lib_fftw3=no)], nfft_lib_fftw3=no)

  # Check if library was found.
  if test "$nfft_lib_fftw3" = "yes"; then
    # Check for user-supplied libdir.
    if test "$with_fftw3" = "yes" -a "$with_fftw3_libdir" = "yes"; then
      # Look for libfftw3.la in /usr/local/lib first
      if test -f "/usr/local/lib/libfftw3.la"; then
        nfft_lib_fftw3_libtool="yes"
      # Look for libfftw3.la in /usr/lib.
      elif test -f "/usr/lib/libfftw3.la"; then
        fftw3_includedir="/usr/include"
        fftw3_includes="-I$fftw3_includedir"
        fftw3_libdir="/usr/lib"
        fftw3_libs="-L$fftw3_libdir"
        nfft_lib_fftw3_libtool="yes"
      else
        nfft_lib_fftw3_libtool="no"
      fi
    else
      # Check user-supplied directory only.
      if test -f "$fftw3_libdir/libfftw3.la"; then
        nfft_lib_fftw3_libtool="yes"
      else
        nfft_lib_fftw3_libtool="no"
      fi
    fi
  fi

  # Restore saved flags.
  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
])
