dnl @synopsis AX_CHECK_DIR(DIR, [ACTION-SUCCESS], [ACTION-FAILURE])
dnl @summary check for directory DIR
dnl @category Misc
dnl
dnl Check whether the directory DIR exists.
dnl
dnl ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
dnl success/failure.
dnl
dnl @version 2008-12-07
dnl @license GPLWithACException
dnl @author Jens Keiner <jens@nfft.org>.
AC_DEFUN([AX_CHECK_DIR],
[
AC_MSG_CHECKING([whether directory $1 exists])
if test -d "$1"; then
  AC_MSG_RESULT(yes)
  m4_default([$2], :)
else
  AC_MSG_RESULT(no)
  m4_default([$3], :)
fi
])dnl AX_CHECK_DIR
