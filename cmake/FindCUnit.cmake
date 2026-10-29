# Copyright (c) 2026 The NFFT Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/>.

# FindCUnit — locate the CUnit unit-testing framework, the NFFT test suite's
# only extra dependency. Probes for the
# <CUnit/CUnit.h> header and a `cunit` library (the suite uses CUnit's
# Automated interface, which needs only -lcunit).
#
# Honors CUnit_ROOT (var or env) and CUnit_INCLUDEDIR / CUnit_LIBDIR.
#
# Provides imported target:
#   CUnit::CUnit
# and sets CUnit_FOUND.

find_path(CUnit_INCLUDE_DIR
  NAMES CUnit/CUnit.h
  HINTS ${CUnit_ROOT} ENV CUnit_ROOT ${CUnit_INCLUDEDIR}
  PATH_SUFFIXES include)

find_library(CUnit_LIBRARY
  NAMES cunit libcunit
  HINTS ${CUnit_ROOT} ENV CUnit_ROOT ${CUnit_LIBDIR}
  PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUnit
  REQUIRED_VARS CUnit_LIBRARY CUnit_INCLUDE_DIR)

if(CUnit_FOUND AND NOT TARGET CUnit::CUnit)
  add_library(CUnit::CUnit UNKNOWN IMPORTED)
  set_target_properties(CUnit::CUnit PROPERTIES
    IMPORTED_LOCATION "${CUnit_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CUnit_INCLUDE_DIR}")
endif()

mark_as_advanced(CUnit_INCLUDE_DIR CUnit_LIBRARY)
