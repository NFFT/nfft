dnl @synopsis AX_CC_MAXOPT
dnl @summary turn on optimization flags for the C compiler
dnl @category C
dnl
dnl Try to turn on "good" C optimization flags for various compilers
dnl and architectures, for some definition of "good".  (In our case,
dnl good for FFTW and hopefully for other scientific codes.  Modify 
dnl as needed.)
dnl
dnl The user can override the flags by setting the CFLAGS environment
dnl variable.  The user can also specify --enable-portable-binary in
dnl order to disable any optimization flags that might result in
dnl a binary that only runs on the host architecture.
dnl
dnl Note also that the flags assume that ANSI C aliasing rules are
dnl followed by the code (e.g. for gcc's -fstrict-aliasing), and that
dnl floating-point computations can be re-ordered as needed.
dnl
dnl Requires macros: AX_CHECK_COMPILER_FLAGS, AX_COMPILER_VENDOR,
dnl                  AX_GCC_ARCHFLAG, AX_GCC_X86_CPUID
dnl
dnl @version 2008-12-02
dnl @license GPLWithACException
dnl @author Steven G. Johnson <stevenj@alum.mit.edu> and Matteo Frigo.
AC_DEFUN([AX_CC_MAXOPT],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_ENABLE(portable-binary, [AC_HELP_STRING([--enable-portable-binary], [disable compiler optimizations that would produce unportable binaries])], 
	acx_maxopt_portable=$enableval, acx_maxopt_portable=no)

# Try to determine "good" native compiler flags if none specified via CFLAGS
if test "$ac_test_CFLAGS" != "set"; then
  CFLAGS=""
  case $ax_cv_c_compiler_vendor in
    dec) CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
      if test "x$acx_maxopt_portable" = xno; then
        CFLAGS="$CFLAGS -arch host"
      fi;;

    sun) CFLAGS="-native -fast -xO5 -dalign"
      if test "x$acx_maxopt_portable" = xyes; then
         CFLAGS="$CFLAGS -xarch=generic"
      fi;;

    hp) CFLAGS="+Oall +Optrs_ansi +DSnative"
      if test "x$acx_maxopt_portable" = xyes; then
        CFLAGS="$CFLAGS +DAportable"
      fi;;

    ibm) if test "x$acx_maxopt_portable" = xno; then
        xlc_opt="-qarch=auto -qtune=auto"
      else
        xlc_opt="-qtune=auto"
      fi
      AX_CHECK_COMPILER_FLAGS($xlc_opt,
        CFLAGS="-O3 -qansialias -w $xlc_opt",
        [CFLAGS="-O3 -qansialias -w"
              echo "******************************************************"
              echo "*  You seem to have the IBM  C compiler.  It is      *"
              echo "*  recommended for best performance that you use:    *"
              echo "*                                                    *"
              echo "*    CFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w *"
              echo "*                      ^^^        ^^^                *"
              echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of *"
              echo "*  CPU you have.  (Set the CFLAGS environment var.   *"
              echo "*  and re-run configure.)  For more info, man cc.    *"
              echo "******************************************************"])
        ;;

    intel) CFLAGS="-O3"
        # Intel seems to have changed the spelling of this flag recently
        icc_ansi_alias="unknown"
	for flag in -ansi-alias -ansi_alias; do
	  AX_CHECK_COMPILER_FLAGS($flag, [icc_ansi_alias=$flag; break])
	done
 	if test "x$icc_ansi_alias" != xunknown; then
            CFLAGS="$CFLAGS $icc_ansi_alias"
        fi
	AX_CHECK_COMPILER_FLAGS(-malign-double, CFLAGS="$CFLAGS -malign-double")
	if test "x$acx_maxopt_portable" = xno; then
	  icc_archflag=unknown
	  icc_flags=""
	  # -xN etcetera are for older versions of icc:
	  case $host_cpu in
	    i686*|x86_64*)
              # icc accepts gcc assembly syntax, so these should work:
	      AX_GCC_X86_CPUID(0)
              AX_GCC_X86_CPUID(1)
	      case $ax_cv_gcc_x86_cpuid_0 in # see AX_GCC_ARCHFLAG
                *:756e6547:*:*) # Intel
                  case $ax_cv_gcc_x86_cpuid_1 in
                    *6a?:*[[234]]:*:*|*6[[789b]]?:*:*:*) icc_flags="-xK";;
                    *f3[[347]]:*:*:*|*f4[1347]:*:*:*) icc_flags="-xP -xN -xW -xK";;
                    *f??:*:*:*) icc_flags="-xN -xW -xK";;
                  esac ;;
              esac ;;
          esac
          # newer icc versions should support -xHost
	  icc_flags="-xHost $icc_flags"
          if test "x$icc_flags" != x; then
            for flag in $icc_flags; do
              AX_CHECK_COMPILER_FLAGS($flag, [icc_archflag=$flag; break])
            done
          fi
          AC_MSG_CHECKING([for icc architecture flag])
	  AC_MSG_RESULT($icc_archflag)
          if test "x$icc_archflag" != xunknown; then
            CFLAGS="$CFLAGS $icc_archflag"
          fi
        fi
	;;
    
    gnu) 
     # default optimization flags for gcc on all systems
     CFLAGS="-O3 -fomit-frame-pointer"

     # -malign-double for x86 systems
     AX_CHECK_COMPILER_FLAGS(-malign-double, CFLAGS="$CFLAGS -malign-double")

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILER_FLAGS(-fstrict-aliasing,
     CFLAGS="$CFLAGS -fstrict-aliasing")

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILER_FLAGS(-ffast-math, CFLAGS="$CFLAGS -ffast-math")

     AX_GCC_ARCHFLAG($acx_maxopt_portable)
     ;;
   apple)
     # default optimization flags for apple on all systems
     AX_CHECK_COMPILER_FLAGS(-O3, CFLAGS="$CFLAGS -O3")
     AX_CHECK_COMPILER_FLAGS(-fomit-frame-pointer, CFLAGS="$CFLAGS -fomit-frame-pointer")
     AX_CHECK_COMPILER_FLAGS(-fPIC, CFLAGS="$CFLAGS -fPIC")

     # -malign-double for x86 systems
     AX_CHECK_COMPILER_FLAGS(-malign-double, CFLAGS="$CFLAGS -malign-double")

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILER_FLAGS(-fstrict-aliasing,CFLAGS="$CFLAGS -fstrict-aliasing")

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILER_FLAGS(-ffast-math, CFLAGS="$CFLAGS -ffast-math")

     AX_GCC_ARCHFLAG($acx_maxopt_portable)
     ;;
  esac

  if test -z "$CFLAGS"; then
	echo ""
	echo "********************************************************"
        echo "* WARNING: Don't know the best CFLAGS for this system  *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
	echo "* (otherwise, a default of CFLAGS=-O3 will be used)    *"
	echo "********************************************************"
	echo ""
        CFLAGS="-O3"
  fi

  AX_CHECK_COMPILER_FLAGS($CFLAGS, [], [
	echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed CFLAGS don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
        echo "********************************************************"
        echo ""
        CFLAGS=""
  ])

fi
])
