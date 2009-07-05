dnl @synopsis AX_APPLE_GCC_ARCHFLAG([PORTABLE],[ACTION-SUCCESS],[ACTION-FAILURE])
dnl @summary find target architecture name for Apple's gcc -arch flag
dnl @category Misc
dnl
dnl This macro tries to guess the "native" arch corresponding to
dnl the target architecture for use with Apple's gcc's -arch flag. If found, the
dnl cache variable $ax_cv_apple_gcc_archflag is set to this flag and 
dnl ACTION-SUCCESS is executed; otherwise $ax_cv_apple_gcc_archflag is
dnl is set to "unknown" and ACTION-FAILURE is executed. The default
dnl ACTION-SUCCESS is to add $ax_cv_apple gcc_archflag to the end of $CFLAGS.
dnl
dnl The user can specify --with-apple-gcc-arch=<arch> in order to override
dnl the macro's choice of architecture, or --without-apple-gcc-arch to
dnl disable this.
dnl
dnl When cross-compiling, or if $CC is not Apple's gcc, then ACTION-FAILURE is
dnl called unless the user specified --with-apple-gcc-arch manually.
dnl
dnl Requires macros: AX_CHECK_COMPILER_FLAGS
dnl
dnl (The main emphasis here is on recent CPUs, on the principle that
dnl  doing high-performance computing on old hardware is uncommon.)
dnl
dnl @version 2008-12-07
dnl @license GPLWithACException
dnl @author Jens Keiner <keiner@math.uni-luebeck.de>.
AC_DEFUN([AX_APPLE_GCC_ARCHFLAG],
[AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_WITH(apple-gcc-arch, [AC_HELP_STRING([--with-apple-gcc-arch=<arch>], 
  [use architecture <arch> for Apple's gcc -arch, instead of guessing])], 
  ax_apple_gcc_arch=$withval, ax_apple_gcc_arch=yes)

AC_MSG_CHECKING([for Apple's gcc architecture flag])
AC_MSG_RESULT([])
AC_CACHE_VAL(ax_cv_apple_gcc_archflag,
[
ax_cv_apple_gcc_archflag="unknown"

if test "$GCC" = yes; then
  if test "x$ax_apple_gcc_arch" = xyes; then
  ax_apple_gcc_arch=""
  if test "$cross_compiling" = no; then
    if test "x[]m4_default([$1],no)" = xyes; then # if we require portable code
      ax_apple_gcc_arch="i386 x86_64 ppc ppc64"
      echo "default arch because of portable code"
    else
      case $host_cpu in
        x86_64*|amd64*|i[[3456]]86*)
          ax_apple_gcc_arch="x86_64 i386"
          ;;
        powerpc*)
           cputype=`((grep cpu /proc/cpuinfo | head -n 1 | cut -d: -f2 | cut -d, -f1 | sed 's/ //g') ; /usr/bin/machine ; /bin/machine; grep CPU /var/run/dmesg.boot | head -n 1 | cut -d" " -f2) 2> /dev/null`
           cputype=`echo $cputype | sed -e 's/ppc//g;s/ *//g'`
           case $cputype in
             *750*|*740[[0-9]]*|*74[[4-5]][[0-9]]*|*74[[0-9]][[0-9]]*) ax_apple_gcc_arch="ppc";;
             *970*|*POWER4*|*power4*|*gq*|*POWER5*|*power5*|*gr*|*gs*) ax_apple_gcc_arch="ppc64";;
             *) ax_apple_gcc_arch="ppc64 ppc";;
           esac
           ;;
        *)
          ax_apple_gcc_arch="x86_64 i386 ppc64 ppc"
          ;;
      esac
    fi # portable code
  fi # not cross-compiling
fi # guess arch

if test "x$ax_apple_gcc_arch" != x -a "x$ax_apple_gcc_arch" != xno; then
  ax_cv_apple_gcc_archflag=""
  for arch in $ax_apple_gcc_arch; do
    AX_CHECK_COMPILER_FLAGS([-arch $arch],[
      saved_CFLAGS="$CFLAGS";
      CFLAGS="$CFLAGS -arch $arch";
      LIBS="$LIBS $fftw3_LIBS"
      AC_MSG_CHECKING([whether linking is possible with -arch $arch]);
      AC_LINK_IFELSE([int main(void){return 0;}],[last_result=yes;AC_MSG_RESULT([yes]);ax_cv_apple_gcc_archflag="$ax_cv_apple_gcc_archflag -arch $arch"],[last_result=no;AC_MSG_RESULT([no])]);
      CFLAGS="$saved_CFLAGS"
    ])
    if test "x$last_result" = "xyes"; then
      break;
    fi
  done
fi

fi # $GCC=yes
])
AC_MSG_CHECKING([for Apple's gcc architecture flag])
AC_MSG_RESULT($ax_cv_apple_gcc_archflag)
if test "x$ax_cv_apple_gcc_archflag" = xunknown; then
  m4_default([$3],:)
else
  m4_default([$2], [CFLAGS="$CFLAGS $ax_cv_apple_gcc_archflag"])
fi
])
