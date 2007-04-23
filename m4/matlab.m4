# $Id: nfsftmex.c 1548 2007-04-18 07:34:23Z keiner $
#
# Copyright (c) 2007 Jens Keiner
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

#serial 2
AC_DEFUN([AC_PROG_MATLAB],
[
  AC_ARG_WITH(matlab,
  [  --with-matlab=DIR    the directory where Matlab is installed ],
  MATLAB_DIR=${withval},MATLAB_DIR=)

  if test -n "${MATLAB_DIR}"; then
    AC_MSG_CHECKING(for Matlab software)
    case $build_os in
      *darwin*)
        MATLAB_CFLAGS="-I${MATLAB_DIR}/extern/include -I${MATLAB_DIR}/simulink/include -fno-common -fexceptions -no-cpp-precomp  -DMATLAB_MEX_FILE -I${MATLAB_DIR}/extern/src";
        MATLAB_LDFLAGS="-bundle -Wl,-flat_namespace -undefined suppress -Wl,-exported_symbols_list,${MATLAB_DIR}/extern/lib/mac/mexFunction.map";
        MATLAB_LDADD="";
        MEXEXT=".mexmac";;
      *linux*)
        MATLAB_CFLAGS="-I${MATLAB_DIR}/extern/include -I${MATLAB_DIR}/simulink/include -fno-common -fexceptions -DMATLAB_MEX_FILE -I${MATLAB_DIR}/extern/src";
        MATLAB_LDFLAGS="-shared -Wl,-flat_namespace -Wl,--version-script,${MATLAB_DIR}/extern/lib/glnx86/mexFunction.map -L${MATLAB_DIR}/bin/glnx86";
        MATLAB_LDADD="-lmx -lmex -lmat -lm";
        MEXEXT=".mexglx";;
      *irix*)
        MATLAB_CFLAGS="-I${MATLAB_DIR}/extern/include -I${MATLAB_DIR}/simulink/include -I${MATLAB_DIR}/extern/src";
        MATLAB_LDFLAGS="-shared -L${MATLAB_DIR}/bin/sgi -L${MATLAB_DIR}/extern/lib/sgi";
        MATLAB_LDADD="-lmx -lmex -lmat -lm";
        MEXEXT=".mexsg";;
      *cygwin*)
        MATLAB_CFLAGS="-DHAVE_WINDOWS -I${MATLAB_DIR}/extern/include -I${MATLAB_DIR}/simulink/include -I${MATLAB_DIR}/extern/src -fno-exceptions -mno-cygwin -DMATLAB_MEX_FILE -DNDEBUG";
        MATLAB_LDFLAGS="-shared -mno-cygwin -W1,--version-script,${MATLAB_DIR}/extern/lib/win32/mexFunction.def -L${MATLAB_DIR}/bin/win32 -W1,--rpath-link,${MATLAB_DIR}/extern/lib/win32,--rpath-link,${MATLAB_DIR}/bin/win32";
        MATLAB_LDADD="-lmx -lmex -lmat -lm";
        MEXEXT=".dll";
        if test ! -e "${MATLAB_DIR}/bin/win32/libmx.a"; then
          cd ${MATLAB_DIR}/bin/win32
          libmx=`dlltool -llibmx.a -d${MATLAB_DIR}/extern/include/libmx.def -Dlibmx.dll`
          cd -
        fi
        if test ! -e "${MATLAB_DIR}/bin/win32/libmex.a"; then
          cd ${MATLAB_DIR}/bin/win32
          libmex=`dlltool -llibmex.a -d${MATLAB_DIR}/extern/include/libmex.def -Dlibmex.dll`
          cd -
        fi
        if test ! -e "${MATLAB_DIR}/bin/win32/libmat.a"; then
          cd ${MATLAB_DIR}/bin/win32
          libmat=`dlltool -llibmat.a -d${MATLAB_DIR}/extern/include/libmat.def -Dlibmat.dll`
          cd -
        fi;;
      *mingw*)
        MATLAB_CFLAGS="-DHAVE_WINDOWS -I${MATLAB_DIR}/extern/include -I${MATLAB_DIR}/simulink/include -fno-exceptions -DMATLAB_MEX_FILE -DNDEBUG";
        MATLAB_LDFLAGS="-shared -W1,--version-script,${MATLAB_DIR}/extern/lib/win32/mexFunction.def -L${MATLAB_DIR}/bin/win32 -W1,--rpath-link,${MATLAB_DIR}/extern/lib/win32,--rpath-link,${MATLAB_DIR}/bin/win32";
        MATLAB_LDADD="-lmx -lmex -lmat -lm";
        MEXEXT=".dll";
        if test ! -e "${MATLAB_DIR}/bin/win32/libmx.a"; then
          cd ${MATLAB_DIR}/bin/win32
          libmx=`dlltool -llibmx.a -d${MATLAB_DIR}/extern/include/libmx.def -Dlibmx.dll`
          cd -
        fi
        if test ! -e "${MATLAB_DIR}/bin/win32/libmex.a"; then
          cd ${MATLAB_DIR}/bin/win32
          libmex=`dlltool -llibmex.a -d${MATLAB_DIR}/extern/include/libmex.def -Dlibmex.dll`
          cd -
        fi
        if test ! -e "${MATLAB_DIR}/bin/win32/libmex.a"; then
          cd ${MATLAB_DIR}/bin/win32
          libmat=`dlltool -llibmat.a -d${MATLAB_DIR}/extern/include/libmat.def -Dlibmat.dll`
          cd -
        fi;;
    esac
    AC_MSG_RESULT([ok])
    AC_SUBST([MATLAB_DIR])
    AC_SUBST([MATLAB_CFLAGS])
    AC_SUBST([MATLAB_LDFLAGS])
    AC_SUBST([MATLAB_LDADD])
    AC_SUBST([MEXEXT])
    AM_CONDITIONAL(HAVE_MATLAB, test "xyes" = "xyes" )
    have_matlab=yes
    AM_CONDITIONAL(HAVE_WINDOWS, test $MEXEXT = ".dll" )
  else
    AM_CONDITIONAL(HAVE_MATLAB, test "xno" = "xyes" )
    AM_CONDITIONAL(HAVE_WINDOWS, test "xno" = "xyes" )
    have_matlab=no
  fi
])
