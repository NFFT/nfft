dnl @synopsis AX_CC_BLOCKS
dnl @summary check if compiler supports blocks
dnl @category Misc
dnl
dnl This macro checks whether the C compiler supports blocks. If it does, the
dnl cache variable $ax_cv_c_cc_blocks is set to "yes", otherwise it is set to 
dnl "no".
dnl
dnl @version 2009-09-14
dnl @license GPLWithACException
dnl @author Jens Keiner <keiner@math.uni-luebeck.de>.
AC_DEFUN([AX_CC_BLOCKS],
[
  AC_CACHE_CHECK([whether $CC supports blocks], ax_cv_c_cc_blocks,
  [
    AC_COMPILE_IFELSE(
    [AC_LANG_SOURCE(
      [AC_LANG_PROGRAM(,[[
void (^my_block)(void);
my_block = ^(void){};
my_block();
      ]])],
      [ax_cv_c_cc_blocks="yes"],
      [ax_cv_c_cc_blocks="no"])
    ])
  ])
])
