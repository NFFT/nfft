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
  AC_REQUIRE([AX_CC_BLOCKS])
  AC_CACHE_VAL(ax_cv_lib_dispatch,
  [
    ax_cv_lib_dispatch="no"
    if test "x$ax_cv_c_cc_blocks" = "xyes"; then
      AC_CHECK_HEADERS([dispatch/dispatch.h],
      [
        AC_MSG_CHECKING([for libdispatch])
        AC_LINK_IFELSE(
          AC_LANG_PROGRAM(
          [[
#include <dispatch/dispatch.h>
          ]],
          [[[
dispatch_queue_t q_default;
q_default = dispatch_get_global_queue(0, 0);
#define COUNT 128
__block double result[COUNT];
dispatch_apply(COUNT, q_default, ^(size_t i){result[i] = 1;});
          ]]]),
        [
          ax_cv_lib_dispatch="yes"
        ],)
      ],)
    fi
    if test "x$ax_cv_lib_dispatch" = "xyes"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
    fi
  ])
])
