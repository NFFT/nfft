# Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# @synopsis AX_PROG_MATLAB
# @summary set up compiler and linker flags for Matlab mex file compilation.
# @category C
#
# @version 2008-12-07
# @license GPLWithACException
# @author Jens Keiner <jens@nfft.org>

AC_DEFUN([AX_PROG_MATLAB],
[
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_REQUIRE([AX_LIB_FFTW3])

  # option to enable mex file compilation for GNU Octave
  AC_ARG_WITH(octave,
    [AS_HELP_STRING([--with-octave=DIR],
      [the directory where GNU Octave is installed])],
    octave_dir=${withval},octave_dir="no")

  # option to enable mex file compilation
  AC_ARG_WITH(matlab,
    [AS_HELP_STRING([--with-matlab=DIR],
      [the directory where Matlab is installed])],
    matlab_dir=${withval},matlab_dir="no")

  AC_ARG_WITH(matlab-arch,
    [AS_HELP_STRING([--with-matlab-arch=DIR],
      [Matlab architecture acronym])],
    matlab_arch=${withval},matlab_arch="yes")

  AC_ARG_ENABLE(matlab-argchecks,
    [AS_HELP_STRING([--enable-matlab-argchecks],
      [Compile Matlab interface with argument checks (recommended) [default=yes]])],
      [ok="$enableval"],
      [ok="yes"])

  if test "x$ok" = "xyes"; then
    AC_DEFINE(MATLAB_ARGCHECKS,1,[Define to enable Matlab argument checks.])
  fi

  AC_ARG_WITH(matlab-fftw3-libdir, [AS_HELP_STRING([--with-matlab-fftw3-libdir=DIR],
  [compile Matlab interface with fftw3 library directory DIR])], matlab_fftw3_lib_dir=$withval, 
    matlab_fftw3_lib_dir="yes")

  AC_ARG_ENABLE(matlab-threads,
    [AS_HELP_STRING([--enable-matlab-threads],
      [Compile Matlab interface with thread support [default same as --enable-openmp]])],
      [matlab_threads="$enableval"],
      [matlab_threads="$enable_threads"])

  if test "x${matlab_dir}" != "xno" -a "x${octave_dir}" != "xno"; then
    AC_MSG_ERROR([The arguments --with-matlab and --with-octave can not be used simultaneously.])
  fi

  AC_MSG_CHECKING([whether to check for Matlab])

  if test "x${matlab_dir}" = "xno"; then
    AC_MSG_RESULT([no])
    ax_prog_matlab="no"
  else
    AC_MSG_RESULT([yes])

    # Matlab root
    AX_CHECK_DIR([${matlab_dir}],[],
      [AC_MSG_ERROR([Please supply a valid path to a Matlab root directory or run configure without the option --with-matlab.])])

    # subdirectories
    matlab_bin_dir="${matlab_dir}/bin"
    matlab_src_dir="${matlab_dir}/extern/src"
    matlab_include_dir="${matlab_dir}/extern/include"
    AX_CHECK_DIR([${matlab_bin_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])
    # More recent versions of Matlab do no longer have the src directory.
    # AX_CHECK_DIR([${matlab_src_dir}],[],
    #   [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])
    AX_CHECK_DIR([${matlab_include_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])

    # architecture and mex file extension
    if test ! "x${matlab_arch}" = "xyes"; then
      AC_MSG_CHECKING([for Matlab architecture])
      AC_MSG_RESULT([${matlab_arch}])
      # mex file extension for architecture
      AC_MSG_CHECKING([for mex file extension])
      case $matlab_arch in
        glnx86) matlab_mexext="mexglx";;
        glnxa64) matlab_mexext="mexa64";;
        mac) matlab_mexext="mexmac";;
        maci) matlab_mexext="mexmaci";;
        maci64) matlab_mexext="mexmaci64";;
        sol64) matlab_mexext="mexs64";;
        win32) matlab_mexext="mexw32";;
        win64) matlab_mexext="mexw64";;
        *) AC_MSG_ERROR([Unsupported or invalid architecture ${matlab_arch}.]);;
      esac
      AC_MSG_RESULT([${matlab_mexext}])
    else 
      matlab_mexext="unknown"
      matlab_arch="unknown"

      # mex file extension, maybe
      for matlab_check_prog_mexext in mexext mexext.sh mexext.bat; do
        unset ac_cv_path_matlab_prog_mexext
        AC_PATH_PROG([matlab_prog_mexext],[$matlab_check_prog_mexext],[no],
          [$matlab_bin_dir$PATH_SEPARATOR$PATH$PATH_SEPARATOR$matlab_bin_dir])
        if test ! "x${matlab_prog_mexext}" = "xno"; then
          AC_MSG_CHECKING([for mex file extension])
	  # Calling mexext.bat from mingw may fail
          if test "x${host_os}" = "xmingw32" -o "x${host_os}" = "xmingw64"; then
            matlab_mexext=unknown
          else
            matlab_mexext=`${matlab_prog_mexext}`
            matlab_mexext=`echo ${matlab_mexext} | tr -d '\r\n'`
          fi  
          AC_MSG_RESULT([${matlab_mexext}])
          break
        fi
      done

      # architecture, maybe
      if test "x${matlab_mexext}" = "xunknown"; then
        # Try guessing the architecture based on host
        case $host in
          *86_64*linux*) matlab_arch_test="glnxa64";;
          *86*linux*) matlab_arch_test="glnx86";;
          *powerpc*darwin*) matlab_arch_test="mac mac64";;
          *86*darwin*) matlab_arch_test="maci maci64";; 
          *solaris*) matlab_arch_test="sol sol64";;
          *cygwin*) matlab_arch_test="win64 win32";;
          *mingw*) matlab_arch_test="win64 win32";;
          *) AC_MSG_ERROR([Cannot guess Matlab architecture based on host type.]);;
        esac
        AC_MSG_CHECKING([for architecture])
        for matlab_arch in $matlab_arch_test; do
          if test -d "${matlab_bin_dir}/${matlab_arch}" -a -f "${matlab_bin_dir}/${matlab_arch}/MATLAB$ac_exeext"; then
            AC_MSG_RESULT([${matlab_arch}])
            break
          fi
          matlab_arch="unkown"
          AC_MSG_RESULT([unknown])
        done
      fi
       
      # mex file extension or architecture found
      if test "x${matlab_mexext}" = "xunknown" -a "x${matlab_arch}" = "xunknown"; then
        AC_MSG_ERROR([Could not determine mex file extension nor Matlab architecture. Please supply a valid architecture flag using the option --with-matlab-arch])
      fi

      if test "x${matlab_arch}" = "xunknown"; then
        # architecture for mex file extension 
        AC_MSG_CHECKING([for architecture])
        case ${matlab_mexext} in
          mexglx) matlab_arch="glnx86";;
          mexa64) matlab_arch="glnxa64";;
          mexmac) matlab_arch="mac";;
          mexmaci) matlab_arch="maci";;
          mexmaci64) matlab_arch="maci64";;
          mexs64) matlab_arch="sol64";;
          mexw32) matlab_arch="win32";;
          mexw64) matlab_arch="win64";;
          *) AC_MSG_ERROR([Unsupported mex file extension ${matlab_mexext}.]);;
        esac
        AC_MSG_RESULT([${matlab_arch}])
      elif test "x${matlab_mexext}" = "xunknown"; then
        AC_MSG_CHECKING([for mex file extension])
        case $matlab_arch in
          glnx86) matlab_mexext="mexglx";;
          glnxa64) matlab_mexext="mexa64";;
          mac) matlab_mexext="mexmac";;
          maci) matlab_mexext="mexmaci";;
          maci64) matlab_mexext="mexmaci64";;
          sol64) matlab_mexext="mexs64";;
          win32) matlab_mexext="mexw32";;
          win64) matlab_mexext="mexw64";;
          *) AC_MSG_ERROR([Unsupported or invalid architecture ${matlab_arch}.]);;
        esac
        AC_MSG_RESULT([${matlab_mexext}])
      fi
    fi

    # add "." to mex file extension
    matlab_mexext=".$matlab_mexext"

    # arch bin dir
    matlab_arch_bin_dir="${matlab_bin_dir}/${matlab_arch}"
    AX_CHECK_DIR([${matlab_arch_bin_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_dir} does not seem to be a valid Matlab root directory.])])

    # dynamic library extension for architecture
    case $matlab_arch in
      glnx86|glnxa64|sol|sol64) matlab_libext=".so";;
      mac|mac64|maci|maci64) matlab_libext=".dylib";;
      win32|win64) matlab_libext=".dll";;
      *) AC_MSG_ERROR([Unsupported or invalid architecture ${matlab_arch}.]);;
    esac

    # libraries
    matlab_LDFLAGS="-L${matlab_arch_bin_dir}"
    saved_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS ${matlab_LDFLAGS}"
    matlab_LIBS=""
    matlab_fftw3_LIBS=""
    matlab_fftw3_LDFLAGS=""
    AC_CHECK_LIB([mx],[mxMalloc],[matlab_LIBS="$matlab_LIBS -lmx"],
      [AC_CHECK_FILE([${matlab_arch_bin_dir}/libmx${matlab_libext}],
      [matlab_LIBS="$matlab_LIBS -lmx"],
      [AC_MSG_ERROR([Needed Matlab library mx not usable.])])],[])
    AC_CHECK_LIB([mex],[mexCallMATLAB],[matlab_LIBS="$matlab_LIBS -lmex"],
      [AC_CHECK_FILE([${matlab_arch_bin_dir}/libmex${matlab_libext}],
      [matlab_LIBS="$matlab_LIBS -lmex"],
      [AC_MSG_ERROR([Needed Matlab library mex not usable.])])],[])
    AC_CHECK_LIB([mat],[matGetVariable],[matlab_LIBS="$matlab_LIBS -lmat"],
      [AC_CHECK_FILE([${matlab_arch_bin_dir}/libmat${matlab_libext}],
      [matlab_LIBS="$matlab_LIBS -lmat"],
      [AC_MSG_ERROR([Needed Matlab library mat not usable.])])],[])

    LDFLAGS="$saved_LDFLAGS"

    matlab_CPPFLAGS="-I${matlab_include_dir} -DMATLAB_DEFAULT_RELEASE=R2017b"

    ax_prog_matlab="yes"

    # Only overwrite Matlab fftw3 lib dir when not explicitly set
    if test "x${matlab_fftw3_lib_dir}" = "xyes"; then
      matlab_fftw3_lib_dir="${matlab_bin_dir}/${matlab_arch}"
    else
      AX_CHECK_DIR([${matlab_fftw3_lib_dir}],[],
      [AC_MSG_ERROR([The directory ${matlab_fftw3_lib_dir} does not seem to be a valid fftw3 library directory.])])
      matlab_fftw3_LDFLAGS="-L$matlab_fftw3_lib_dir"
    fi

    saved_LIBS="$LIBS"
    saved_LDFLAGS="$LDFLAGS"

    for matlab_fftw3_lib_name in mwfftw3${PREC_SUFFIX} :libmwfftw3${PREC_SUFFIX}.so.3 fftw3${PREC_SUFFIX}; do
      matlab_fftw3_LIBS="-l${matlab_fftw3_lib_name}"
      LIBS="-l${matlab_fftw3_lib_name} $LIBS"
      LDFLAGS="-L$matlab_fftw3_lib_dir $LDFLAGS"
      AC_MSG_CHECKING([for Matlab fftw3 library])
      AC_LINK_IFELSE([AC_LANG_CALL([], [fftw${PREC_SUFFIX}_execute])], [ax_matlab_lib_fftw3=yes],[ax_matlab_lib_fftw3=no])
      AC_MSG_RESULT([$ax_matlab_lib_fftw3])

      if test "x$ax_matlab_lib_fftw3" = "xno"; then
        matlab_fftw3_LIBS="-l${matlab_fftw3_lib_name} -lm"
        LIBS="$matlab_fftw3_LIBS $saved_LIBS"
        AC_MSG_CHECKING([for Matlab fftw3 library (-lm set)])
        AC_LINK_IFELSE([AC_LANG_CALL([], [fftw${PREC_SUFFIX}_execute])], [ax_matlab_lib_fftw3=yes],[ax_matlab_lib_fftw3=no])
        AC_MSG_RESULT([$ax_matlab_lib_fftw3])
      fi

      if test "x$ax_matlab_lib_fftw3" = "xno"; then
        matlab_fftw3_LIBS="-l${matlab_fftw3_lib_name} -pthread -lm"
        LIBS="$matlab_fftw3_LIBS $saved_LIBS"
        AC_MSG_CHECKING([for Matlab fftw3 library (-lpthread -lm set)])
        AC_LINK_IFELSE([AC_LANG_CALL([], [fftw${PREC_SUFFIX}_execute])], [ax_matlab_lib_fftw3=yes],[ax_matlab_lib_fftw3=no])
        AC_MSG_RESULT([$ax_matlab_lib_fftw3])
      fi

      if test "x$ax_matlab_lib_fftw3" = "xno"; then
        continue
      fi

      if test "x$matlab_threads" = "xyes"; then
        LIBS="$matlab_fftw3_LIBS $saved_LIBS"
        AC_MSG_CHECKING([for Matlab combined fftw3 library with thread support])
        AC_LINK_IFELSE([AC_LANG_CALL([], [fftw_init_threads])], [ax_matlab_lib_fftw3_threads=yes],[ax_matlab_lib_fftw3_threads=no])
        AC_MSG_RESULT([$ax_matlab_lib_fftw3_threads])

        if test "x$ax_matlab_lib_fftw3_threads" = "xno"; then
          ax_matlab_lib_fftw3_threads="yes"
          LIBS="$matlab_fftw3_LIBS -lpthread -lm $saved_LIBS"
          AC_MSG_CHECKING([for Matlab combined fftw3 library with thread support (-lpthread -lm set)])
          AC_LINK_IFELSE([AC_LANG_CALL([], [fftw${PREC_SUFFIX}_init_threads])], [matlab_fftw3_LIBS="$matlab_fftw3_LIBS -lpthread -lm"],[ax_matlab_lib_fftw3_threads=no])
          AC_MSG_RESULT([$ax_matlab_lib_fftw3_threads])
        fi

        if test "x$ax_matlab_lib_fftw3_threads" = "xno"; then
          ax_matlab_lib_fftw3_threads="yes"
          LIBS="-l${matlab_fftw3_lib_name}_threads $matlab_fftw3_LIBS $saved_LIBS"
          AC_MSG_CHECKING([for Matlab fftw3 library with thread support])
          AC_LINK_IFELSE([AC_LANG_CALL([], [fftw${PREC_SUFFIX}_init_threads])], [matlab_fftw3_LIBS="-l${matlab_fftw3_lib_name}_threads $matlab_fftw3_LIBS"],[ax_matlab_lib_fftw3_threads=no])
          AC_MSG_RESULT([$ax_matlab_lib_fftw3_threads])
        fi

        if test "x$ax_matlab_lib_fftw3_threads" = "xno"; then
          ax_matlab_lib_fftw3_threads="yes"
          LIBS="-l${matlab_fftw3_lib_name}_threads -lpthread $matlab_fftw3_LIBS $saved_LIBS"
          AC_MSG_CHECKING([for Matlab fftw3 library with thread support (-lpthread set)])
          AC_LINK_IFELSE([AC_LANG_CALL([], [fftw${PREC_SUFFIX}_init_threads])], [matlab_fftw3_LIBS="-l${matlab_fftw3_lib_name}_threads -lpthread $matlab_fftw3_LIBS"],[ax_matlab_lib_fftw3_threads=no])
          AC_MSG_RESULT([$ax_matlab_lib_fftw3_threads])
        fi

        if test "x$ax_matlab_lib_fftw3_threads" = "xno"; then
          continue
        fi
      fi

      if test "x$ax_matlab_lib_fftw3" != "xno"; then
        break
      fi
    done

    if test "x$ax_matlab_lib_fftw3" = "xno"; then
      AC_MSG_ERROR([You don't seem to have installed the FFTW 3 libray for the NFFT Matlab interface.])
    fi

    if test "x$matlab_threads" = "xyes" -a "x$ax_matlab_lib_fftw3_threads" = "xno"; then
      AC_MSG_ERROR([You don't seem to have installed the threaded FFTW 3 libray for the NFFT Matlab interface.])
    fi

    LIBS="$saved_LIBS"
    LDFLAGS="$saved_LDFLAGS"
  fi


  AC_MSG_CHECKING([whether to check for GNU Octave])

  if test "x${octave_dir}" = "xno"; then
    AC_MSG_RESULT([no])
    ax_prog_octave="no"
  else
    AC_MSG_RESULT([yes])
    ax_prog_octave="yes"

    if test "x${octave_dir}" = "xyes"; then
      octave_dir=""
    fi

    # Modified code from RcppOctave
    AC_MSG_CHECKING([Octave custom binary path specification])
    if test [ -n "$octave_dir" ] ; then # passed as an option
      AC_MSG_RESULT([$octave_dir [[from configure option --with-octave]]])
      if test [ -n "${OCTAVE_PATH}" ] ; then
	AC_MSG_NOTICE([overriding environment variable \$OCTAVE_PATH])
      fi
      OCTAVE_PATH="$octave_dir"
    elif test [ -n "${OCTAVE_PATH}" ] ; then
      AC_MSG_RESULT([${OCTAVE_PATH} [[from environment variable OCTAVE_PATH]]])
    else
      AC_MSG_RESULT([none])
    fi

    # build lookup path for octave-config
    AS_IF([test -n "${OCTAVE_PATH}"], [
	if test [ -f "$OCTAVE_PATH" ] ; then # path is a file: use parent directory
	  OCTAVE_PATH=`AS_DIRNAME(["$OCTAVE_PATH"])`
	fi
	OCTAVE_LOOKUP_PATH="${OCTAVE_PATH}${PATH_SEPARATOR}${OCTAVE_PATH}/bin${PATH_SEPARATOR}${OCTAVE_PATH}/mingw64/bin${PATH_SEPARATOR}${OCTAVE_PATH}/mingw32/bin"
	],[	
	OCTAVE_LOOKUP_PATH="$PATH"
	AC_MSG_NOTICE([using Octave binary path from \$PATH])
	]
    )

    AC_PATH_PROG([OCTAVE_CONFIG], [octave-config], [], [${OCTAVE_LOOKUP_PATH}])
    AC_PATH_PROG([OCTAVE_MKOCTFILE], [mkoctfile], [], [${OCTAVE_LOOKUP_PATH}])
    AC_PATH_PROG([OCTAVE_CLI], [octave-cli], [], [${OCTAVE_LOOKUP_PATH}])

    if test "${OCTAVE_CLI}" != ""; then
      octave_cli="${OCTAVE_CLI}"
    else
      octave_cli="${OCTAVE_PATH}"/bin/octave-cli
    fi

    AC_ARG_WITH(octave-libdir, [AS_HELP_STRING([--with-octave-libdir=DIR],
    [compile with Octave library directory DIR])], octave_lib_dir=$withval, 
      octave_lib_dir="yes")

    AC_ARG_WITH(octave-includedir, [AS_HELP_STRING([--with-octave-includedir=DIR],
    [compile with octave include directory DIR])], octave_include_dir=$withval, 
      octave_include_dir="yes")

    if test "x${octave_include_dir}" = "xyes"; then
      if test "${OCTAVE_CONFIG}" != ""; then	
        AC_MSG_CHECKING([Octave includes directory])
        octave_include_dir=`${OCTAVE_CONFIG} --print OCTINCLUDEDIR`
        if test "x${host_os}" = "xmingw32" -o "x${host_os}" = "xmingw64"; then
          octave_include_dir=`cygpath -u ${octave_include_dir}`
        fi
	AC_MSG_RESULT([${octave_include_dir}])
      else
        octave_include_dir=""
      fi
    fi

    if test "x${octave_lib_dir}" = "xyes"; then 
      if test "${OCTAVE_CONFIG}" != ""; then	
        AC_MSG_CHECKING([Octave libraries directory])
	octave_lib_dir=`${OCTAVE_CONFIG} --print OCTLIBDIR`
        if test "x${host_os}" = "xmingw32" -o "x${host_os}" = "xmingw64"; then
	  octave_lib_dir=`cygpath -u ${octave_lib_dir}`
        fi
	AC_MSG_RESULT([${octave_lib_dir}])
      else
        octave_lib_dir=""
      fi
    fi

    if test [ -n "$octave_include_dir" ]; then
      matlab_CPPFLAGS="-I${octave_include_dir}"
    else
      matlab_CPPFLAGS=""
    fi

    if test [ -n "$octave_lib_dir" ]; then
      matlab_LDFLAGS="-L${octave_lib_dir}"
    else
      matlab_LDFLAGS=""
    fi

    matlab_fftw3_LDFLAGS="$fftw3_LDFLAGS"
    if test "x$matlab_threads" = "xyes"; then
      matlab_fftw3_LIBS="$fftw3_threads_LIBS"
    else
      matlab_fftw3_LIBS="$fftw3_LIBS"
    fi
    matlab_mexext=".mex"

    saved_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$matlab_CPPFLAGS"
    AC_CHECK_HEADER([mex.h], [], [AC_MSG_ERROR([Needed mex.h not found.])])
    CPPFLAGS="$saved_CPPFLAGS"

    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    if test "${OCTAVE_MKOCTFILE}" != ""; then	
      AC_MSG_CHECKING([Octave liboctave flag])
      octave_liboctave=`${OCTAVE_MKOCTFILE} --print LIBOCTAVE`
      AC_MSG_RESULT([${octave_liboctave}])

      AC_MSG_CHECKING([Octave libinterp flag])
      octave_liboctinterp=`${OCTAVE_MKOCTFILE} --print LIBOCTINTERP`
      AC_MSG_RESULT([${octave_liboctinterp}])

      matlab_LIBS=""

      if test [ -n "${octave_liboctave}"]; then
        LDFLAGS="${saved_LDFLAGS} ${matlab_LDFLAGS}"
        LIBS="${saved_LIBS} ${octave_liboctave}"
        AC_MSG_CHECKING([for usable ${octave_liboctave}])
        AC_LINK_IFELSE([AC_LANG_CALL([], [octave_handle_signal])], [
          AC_MSG_RESULT([yes])
          matlab_LIBS="${octave_liboctave}"
          ],[AC_MSG_ERROR([no])])
      fi

      if test [ -n "${octave_liboctinterp}"]; then
        LDFLAGS="${saved_LDFLAGS} ${matlab_LDFLAGS}"
        LIBS="${saved_LIBS} ${octave_liboctinterp} ${matlab_LIBS}"
        AC_MSG_CHECKING([for usable ${octave_liboctinterp}])
        AC_LINK_IFELSE([AC_LANG_CALL([], [mexCallMATLAB])], [
          AC_MSG_RESULT([yes])
          matlab_LIBS="${octave_liboctinterp} ${matlab_LIBS}"
          ],[
          AC_MSG_ERROR([no])
          ])
      elif test [ -n "${octave_liboctave}"]; then
        LDFLAGS="$saved_LDFLAGS ${matlab_LDFLAGS}"
        LIBS="${saved_LIBS} ${matlab_LIBS}"
        AC_MSG_CHECKING([for usable ${octave_liboctave}])
        AC_LINK_IFELSE([AC_LANG_CALL([], [mexCallMATLAB])], [
          AC_MSG_RESULT([yes])
          ],[
          AC_MSG_ERROR([no])
          ])
      fi
    else
      matlab_LIBS="-loctinterp -loctave"
    fi

    LDFLAGS="$saved_LDFLAGS ${matlab_LDFLAGS}"
    LIBS="${saved_LIBS} ${matlab_LIBS}"
    AC_MSG_CHECKING([for usable Octave MEX interface])
    AC_LINK_IFELSE([AC_LANG_CALL([], [mexCallMATLAB])], [AC_MSG_RESULT([yes])],[AC_MSG_ERROR([no])])

    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

  fi


  AM_CONDITIONAL(HAVE_OCTAVE, test "x$ax_prog_octave" = "xyes" )
  AM_CONDITIONAL(HAVE_MATLAB, test "x$ax_prog_matlab" = "xyes" -o "x$ax_prog_octave" = "xyes" )
  AM_CONDITIONAL(HAVE_MATLAB_THREADS, test "x$matlab_threads" = "xyes")
  AC_SUBST(matlab_CPPFLAGS)
  AC_SUBST(matlab_LIBS)
  AC_SUBST(matlab_LDFLAGS)
  AC_SUBST(matlab_bin_dir)
  AC_SUBST(matlab_mexext)
  AC_SUBST(matlab_fftw3_LIBS)
  AC_SUBST(matlab_fftw3_LDFLAGS)
  AC_SUBST(octave_dir)
  AC_SUBST(octave_cli)
])
