AC_DEFUN([AX_NFFT_MODULE],
[
  AC_ARG_ENABLE($1, [AC_HELP_STRING([--enable-]$1,[build ]$2[ module (]$3[)])], 
    ok=$enableval, ok=$nfft_module_default)
  if m4_default($6,test "x$ok" = "xyes"); then
    AC_DEFINE(HAVE_$2, 1, [Define to enable ]$2[module.])
    HAVE_$2="#define HAVE_$2 1"
    $4
  else
    HAVE_$2="#undef HAVE_$2"
    $5
  fi
  AM_CONDITIONAL(HAVE_$2, m4_default($6,test "x$ok" = "xyes"))
  AC_SUBST(HAVE_$2)
])