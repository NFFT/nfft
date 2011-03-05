dnl @synopsis AX_LIB_DISPATCH
dnl @summary check if libdispatch is available
dnl @category Misc
dnl
dnl This macro checks whether the C library libdispatch is available.
dnl
dnl @version 2009-09-14
dnl @license GPLWithACException
dnl @author Jens Keiner <keiner@math.uni-luebeck.de>.
AC_DEFUN([AX_LIB_DISPATCH],
[
  AC_ARG_ENABLE(libdispatch, [AC_HELP_STRING([--enable-libdispatch], 
    [enable support for libdispatch (EXPERIMENTAL AND UNSUPPORTED)])], want_dispatch=$enableval, want_dispatch=no)
  if test "x$want_dispatch" = "xyes"; then
    AC_CACHE_VAL(ax_cv_lib_dispatch,
    [
      ax_cv_lib_dispatch="no"
      AX_CC_BLOCKS
      if test "x$ax_cv_c_cc_blocks" = "xyes"; then
        AC_CHECK_HEADERS([dispatch/dispatch.h],
        [
          AC_LINK_IFELSE(
            [AC_LANG_SOURCE(
              [AC_LANG_PROGRAM(
              [[
#include <dispatch/dispatch.h>
              ]],
              [[[
  dispatch_queue_t q_default;
  q_default = dispatch_get_global_queue(0, 0);
  #define COUNT 128
  __block double result[COUNT];
  dispatch_apply(COUNT, q_default, ^(size_t i){result[i] = 1;});
              ]]])])],
            [
            ax_cv_lib_dispatch="yes"
            ],[ax_lib_dispatch_reason="could not link against libdispatch"])
          ],[ax_lib_dispatch_reason="dispatch/dispatch.h not found or unusable"])
      else
        ax_lib_dispatch_reason="$CC does not support blocks"
      fi
    ])
    AC_MSG_CHECKING([for libdispatch])
    if test "x$ax_cv_lib_dispatch" = "xyes" -a "x$want_dispatch" = "xyes"; then
      AC_MSG_RESULT([yes])
      AC_DEFINE(HAVE_LIBDISPATCH,1,[Define to enable code that uses libdispatch.])
    else
      AC_MSG_RESULT([no ($ax_lib_dispatch_reason)])
    fi
  fi
])

