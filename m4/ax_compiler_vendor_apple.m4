dnl @synopsis AX_COMPILER_VENDOR_APPLE
dnl @summary check if apple is the compiler vendor of the C compiler
dnl @category C
dnl
dnl Determine if Apple is the vendor of the C compiler. The result is
dnl returned in the variable $ax_c_compiler_vendor_apple.
dnl
dnl @version 2008-11-30
dnl @license GPLWithACException
dnl @author Jens Keiner <keiner@math.uni-luebeck.de>

AC_DEFUN([AX_COMPILER_VENDOR_APPLE],
[
AC_MSG_CHECKING([whether C compiler is Apple's gcc])
AC_RUN_IFELSE(
  [
    AC_LANG_PROGRAM(
      [],
      [
        #if defined(__APPLE_CC__) 
          return 0;
        #else 
          return 1;
        #endif
      ])
  ],
  [
    ax_c_compiler_vendor_apple="yes"
    AC_MSG_RESULT([yes])
  ],
  [
    ax_c_compiler_vendor_apple="no"
    AC_MSG_RESULT([no])
  ])
])

