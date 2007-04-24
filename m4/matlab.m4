dnl $Id$
dnl
dnl Copyright (c) 2007 Jens Keiner
dnl
dnl This program is free software; you can redistribute it and/or modify it 
dnl under the terms of the GNU General Public License as published by the Free 
dnl Software Foundation; either version 2 of the License, or (at your option) 
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful, but WITHOUT
dnl ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
dnl FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
dnl more details.
dnl
dnl You should have received a copy of the GNU General Public License along with
dnl this program; if not, write to the Free Software Foundation, Inc., 51
dnl Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
dnl
dnl @synopsis AX_PROG_MATLAB
dnl @summary set up compiler and linker flags for Matlab mex file compilation.
dnl @category C
dnl
dnl @version 2007-04-24
dnl @license GPLWithACException
dnl @author Jens Keiner

AC_DEFUN([AX_PROG_MATLAB],
[
  AC_REQUIRE([AC_CANONICAL_HOST])
  
  dnl Add configure option to enable mex file compilation.
  AC_ARG_WITH(matlab,
  [  --with-matlab=DIR    the directory where Matlab is installed ],
  MATLAB_DIR=${withval},MATLAB_DIR=)

  dnl Test if Matlab directory existent.
  if test -n "${MATLAB_DIR}"; then
    AC_MSG_CHECKING(for Matlab)
    missing=
    MATLAB_CFLAGS_COMMON="-I${MATLAB_DIR}/extern/include -I${MATLAB_DIR}/extern/src -DMATLAB_MEX_FILE"
    MATLAB_LDADD_COMMON="-lmx -lmex -lmat -lm"
    case $host in
      powerpc-*darwin*)
        dnl Mac (PowerPC)
        AX_CHECK_FILE([${MATLAB_DIR}/bin/mac/libmat.dylib])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/mac/libmex.dylib])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/mac/libmx.dylib])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mat.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/matrix.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mex.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/src/mexversion.c])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/lib/mac/mexFunction.map])
        AC_DEFINE([HAVE_MEXVERSION_C],[1],[Define to have the file mexversion.c])
        MATLAB_CFLAGS="$MATLAB_CFLAGS_COMMON -fno-common -fexceptions -no-cpp-precomp";
        MATLAB_LDFLAGS="-bundle -Wl,-flat_namespace -undefined suppress -Wl,-exported_symbols_list,${MATLAB_DIR}/extern/lib/mac/mexFunction.map";
        MATLAB_LDADD="";
        MEXEXT=".mexmac";;
      i686-*darwin*)
        dnl Mac (Intel)
        dnl FIXME Needs to be tested
        AX_CHECK_FILE([${MATLAB_DIR}/bin/maci/libmat.dylib])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/maci/libmex.dylib])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/maci/libmx.dylib])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mat.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/matrix.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mex.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/src/mexversion.c])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/lib/maci/mexFunction.map])
        AC_DEFINE([HAVE_MEXVERSION_C],[1],[Define to have the file mexversion.c])
        MATLAB_CFLAGS="$MATLAB_CFLAGS_COMMON -fno-common -fexceptions -no-cpp-precomp";
        MATLAB_LDFLAGS="-bundle -Wl,-flat_namespace -undefined suppress -Wl,-exported_symbols_list,${MATLAB_DIR}/extern/lib/mac/mexFunction.map";
        MATLAB_LDADD="";
        MEXEXT=".mexmaci";;
      *86-*linux*)   
        dnl Linux (x86, 32 bit)    
        dnl FIXME Add 64 bit Linux target
        AX_CHECK_FILE([${MATLAB_DIR}/bin/glnx86/libmat.so])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/glnx86/libmex.so])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/glnx86/libmx.so])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mat.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/matrix.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mex.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/src/mexversion.c])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/lib/glnx86/mexFunction.map])               
        AC_DEFINE([HAVE_MEXVERSION_C],[1],[Define to have the file mexversion.c])
        MATLAB_CFLAGS="$MATLAB_CFLAGS_COMMON -fno-common -fexceptions";
        MATLAB_LDFLAGS="-shared -Wl,-flat_namespace -Wl,--version-script,${MATLAB_DIR}/extern/lib/glnx86/mexFunction.map -L${MATLAB_DIR}/bin/glnx86";
        MATLAB_LDADD="$MATLAB_LDFLAGS_COMMON";
        MEXEXT=".mexglx";;
      *irix*)
        dnl IRIX (mips)
        dnl FIXME Needs to be tested
        AX_CHECK_FILE([${MATLAB_DIR}/bin/sgi/libmat.so])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/sgi/libmex.so])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/sgi/libmx.so])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mat.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/matrix.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mex.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/src/mexversion.c])
        AC_DEFINE([HAVE_MEXVERSION_C],[1],[Define to have the file mexversion.c])
        MATLAB_CFLAGS="$MATLAB_CFLAGS_COMMON";
        MATLAB_LDFLAGS="-shared -L${MATLAB_DIR}/bin/sgi -L${MATLAB_DIR}/extern/lib/sgi";
        MATLAB_LDADD="$MATLAB_LDFLAGS_COMMON";
        MEXEXT=".mexsg";;
      *solaris*)
        dnl FIXME Matlab on Solaris only supported on 64bit systems       
        dnl FIXME Needs to be tested
        AX_CHECK_FILE([${MATLAB_DIR}/bin/sol/libmat.so])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/sol/libmex.so])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/sol/libmx.so])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mat.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/matrix.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mex.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/src/mexversion.c])
        AC_DEFINE([HAVE_MEXVERSION_C],[1],[Define to have the file mexversion.c])
        MATLAB_CFLAGS="$MATLAB_CFLAGS_COMMON";
        MATLAB_LDFLAGS="-shared -L${MATLAB_DIR}/bin/sol -L${MATLAB_DIR}/extern/lib/sol";
        MATLAB_LDADD="$MATLAB_LDFLAGS_COMMON";
        MEXEXT=".mexsg";;
      *cygwin*)
        dnl Windows (Cygwin)
        dnl FIXME Needs to be tested
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
        fi
        AX_CHECK_FILE([${MATLAB_DIR}/bin/win32/libmat.a])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/win32/libmex.a])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/win32/libmx.a])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mat.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/matrix.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mex.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/lib/win32/mexFunction.map])               
        MATLAB_CFLAGS="$MATLAB_CFLAGS_COMMON -fno-exceptions -mno-cygwin -DNDEBUG -DHAVE_WINDOWS";
        MATLAB_LDFLAGS="-shared -mno-cygwin -W1,--version-script,${MATLAB_DIR}/extern/lib/win32/mexFunction.def -L${MATLAB_DIR}/bin/win32 -W1,--rpath-link,${MATLAB_DIR}/extern/lib/win32,--rpath-link,${MATLAB_DIR}/bin/win32";
        MATLAB_LDADD="$MATLAB_LDFLAGS_COMMON";
        MEXEXT=".dll";;
      *mingw*)
        dnl Windows (MinGW)
        dnl FIXME Needs to be tested
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
        fi
        AX_CHECK_FILE([${MATLAB_DIR}/bin/win32/libmat.a])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/win32/libmex.a])
        AX_CHECK_FILE([${MATLAB_DIR}/bin/win32/libmx.a])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mat.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/matrix.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/include/mex.h])
        AX_CHECK_FILE([${MATLAB_DIR}/extern/lib/win32/mexFunction.map])               
        MATLAB_CFLAGS="$MATLAB_CFLAGS_COMMON -fno-exceptions -DNDEBUG -DHAVE_WINDOWS";
        MATLAB_LDFLAGS="-shared -W1,--version-script,${MATLAB_DIR}/extern/lib/win32/mexFunction.def -L${MATLAB_DIR}/bin/win32 -W1,--rpath-link,${MATLAB_DIR}/extern/lib/win32,--rpath-link,${MATLAB_DIR}/bin/win32";
        MATLAB_LDADD="$MATLAB_LDFLAGS_COMMON";
        MEXEXT=".dll";;
    esac
    AC_SUBST([MATLAB_DIR])
    AC_SUBST([MATLAB_CFLAGS])
    AC_SUBST([MATLAB_LDFLAGS])
    AC_SUBST([MATLAB_LDADD])
    AC_SUBST([MEXEXT])
    AM_CONDITIONAL(HAVE_MATLAB, test "xyes" = "xyes" )
    AM_CONDITIONAL(MEX_DLL_HACK, test $MEXEXT = ".dll" )
    AX_PROG_MATLAB_RESULT()
  else
    AC_MSG_CHECKING(for Matlab)
    AM_CONDITIONAL(HAVE_MATLAB, test "xno" = "xyes" )
    AM_CONDITIONAL(MEX_DLL_HACK, test "xno" = "xyes" )
    AC_MSG_RESULT([mex file compilation is turned off])
  fi
])

AC_DEFUN([AX_CHECK_FILE],
[
  if ! test -e $1; then missing="$missing\n$1"; fi
])

AC_DEFUN([AX_PROG_MATLAB_RESULT],
[
  if test -n "$missing"; then
    AC_MSG_RESULT([no (missing files follow)])
    echo $missing;
  else
    AC_MSG_RESULT([ok])
  fi                     
])
