dnl
dnl  NFFT_LIB_FFTW3
dnl
dnl  Check configure options and assign variables related to the fftw3 library.
dnl
dnl  If we find the library, set the shell variable `nfft_lib_fftw3' to `yes'.
dnl  Otherwise, set `nfft_lib_fftw3' to `no'.

AC_DEFUN([NFFT_LIB_FFTW3],
[
  AC_ARG_WITH(fftw3, [AC_HELP_STRING([--with-fftw3=ARG],
  [compile with fftw3 in PATH])], with_fftw3=$withval, with_fftw3="yes")

  AC_MSG_CHECKING([whether to check for fftw3 library])

  if test "${with_fftw3}" = "no"; then
    AC_MSG_RESULT([no])
    nfft_lib_fftw3=no
  else
    AC_MSG_RESULT([yes])
    saved_LDFLAGS="$LDFLAGS"
    saved_CFLAGS="$CFLAGS"

    dnl If the user doesn't specify a (valid) directory (or he doesn't supply a
    dnl --with-fftw3 option at all), we want to look in the default directories:
    dnl /usr and /usr/local. However, the compiler always looks in
    dnl /usr/{lib,include} anyway, so we only need to look in /usr/local

    if test ! -d ${with_fftw3}; then
      AC_MSG_NOTICE([Looking in default locations])
      with_fftw3="/usr/local"
    fi

    FFTW3_INCLUDES="-I${with_fftw3}/include"
    FFTW3_LIBDIR="${with_fftw3}/lib"
    FFTW3_LIBS="-L$FFTW3_LIBDIR"
    CFLAGS="$CFLAGS $FFTW3_INCLUDES"
    LDFLAGS="$LDFLAGS $FFTW3_LIBS"

    AC_CHECK_HEADER(fftw3.h,
      [AC_CHECK_LIB(fftw3, fftw_execute,
        nfft_lib_fftw3=yes,
        nfft_lib_fftw3=no)], nfft_lib_fftw3=no)

    if test ! -f "$FFTW3_LIBDIR/libfftw3.la"; then
      AC_MSG_NOTICE([fftw3 does not seem to be installed as a libtool library])
      nfft_lib_fftw3_libtool="no"
    else
      nfft_lib_fftw3_libtool="yes"
    fi

    CFLAGS="$saved_CFLAGS"
    LDFLAGS="$saved_LDFLAGS"
  fi
])
