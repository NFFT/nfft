#serial 3

dnl From Jens Keiner

dnl Define HAVE_ZERO_BYTE_DOUBLE if `double' with value zero consits of bytes containing zeros

AC_DEFUN([AC_ZERO_BYTES_DOUBLE],
[
  AC_MSG_CHECKING(if doubles of value zero consist of bytes containing zeros)
  AC_CACHE_VAL(ac_cv_have_zero_byte_double,
    [AC_TRY_RUN(
	  [#include "stdlib.h"
    #include "string.h"
    int main(void)
    {
      double x;
      memset(&x,0,sizeof(double));
      if (x == 0.0)
        exit(0);
      else  
        exit(1);
    }
	  ],
	  ac_cv_have_zero_byte_double=yes, 
	  ac_cv_have_zero_byte_double=no, 
	  ac_cv_have_zero_byte_double=cross)
	])	 
  AC_MSG_RESULT([$ac_cv_have_zero_byte_double])
if test $ac_cv_have_zero_byte_double = yes; then
  AC_DEFINE(HAVE_ZERO_BYTES_DOUBLE, 1, "")
fi])


dnl From Jim Meyering

dnl Define HAVE_STRUCT_UTIMBUF if `struct utimbuf' is declared --
dnl usually in <utime.h>.
dnl Some systems have utime.h but don't declare the struct anywhere.

AC_DEFUN(jm_CHECK_TYPE_STRUCT_UTIMBUF,
[
  AC_CHECK_HEADERS(utime.h)
  AC_REQUIRE([AC_HEADER_TIME])
  AC_CACHE_CHECK([for struct utimbuf], fu_cv_sys_struct_utimbuf,
    [AC_TRY_COMPILE(
      [
#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif
#ifdef HAVE_UTIME_H
# include <utime.h>
#endif
      ],
      [static struct utimbuf x; x.actime = x.modtime;],
      fu_cv_sys_struct_utimbuf=yes,
      fu_cv_sys_struct_utimbuf=no)
    ])

  if test $fu_cv_sys_struct_utimbuf = yes; then
    AC_DEFINE_UNQUOTED(HAVE_STRUCT_UTIMBUF, 1,
[Define if struct utimbuf is declared -- usually in <utime.h>.
   Some systems have utime.h but don't declare the struct anywhere. ])
  fi
])
