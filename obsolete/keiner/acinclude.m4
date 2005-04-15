# From Jens Keiner
# Define HAVE_ZERO_BYTE_DOUBLE if `double' with value zero consits of bytes containing zeros

AC_DEFUN([AC_ZERO_BYTES_DOUBLE],
[
  AC_MSG_CHECKING(if doubles of value zero consist of bytes containing zeros)
  AC_CACHE_VAL(ac_cv_have_zero_byte_double,
    [AC_TRY_RUN(
	  [#include <stdlib.h>
    #include <string.h>
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
