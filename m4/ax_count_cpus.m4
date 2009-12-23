# ===========================================================================
#          http://www.nongnu.org/autoconf-archive/ax_count_cpus.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_COUNT_CPUS
#
# DESCRIPTION
#
#   Attempt to count the number of processors present on the machine. If the
#   detection fails, then a value of 1 is assumed.
#
#   The value is placed in the CPU_COUNT variable.
#
# LICENSE
#
#   Copyright (c) 2008 Michael Paul Bailey <jinxidoru@byu.net>
#   Copyright (c) 2008 Christophe Tournayre <turn3r@users.sourceforge.net>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_COUNT_CPUS], [
    AC_REQUIRE([AC_PROG_EGREP])
    AC_MSG_CHECKING(the number of available CPUs)
    cpu_count="0"

    #On MacOS
    if test -x /usr/sbin/sysctl -a `/sbin/sysctl -a 2>/dev/null| grep -c hw.cpu`; then
        cpu_count=`/usr/sbin/sysctl -n hw.ncpu`
    fi

    #On Linux
    if test "x$CPU_COUNT" = "x0" -a -e /proc/cpuinfo; then
        cpu_count=`$EGREP -c '^processor' /proc/cpuinfo`
    fi

    if test "x$cpu_count" = "x0"; then
        cpu_count="1"
        AC_MSG_RESULT( [unable to detect (assuming 1)] )
    else
        AC_MSG_RESULT( $cpu_count )
    fi
])
