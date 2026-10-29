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

# Julia interface (double + shared only).
nfft_add_julia_module(nfft libnfftjulia.c)
if(HAVE_NFCT)
  nfft_add_julia_module(nfct libnfctjulia.c)
endif()
if(HAVE_NFST)
  nfft_add_julia_module(nfst libnfstjulia.c)
endif()
if(HAVE_NFSFT)
  nfft_add_julia_module(nfsft libnfsftjulia.c)
endif()
if(NFFT_ENABLE_APPLICATIONS)
  if(NFFT_ENABLE_OPENMP)
    nfft_add_julia_module(fastsum libfastsumjulia.c
      EXTRA_INC ${PROJECT_SOURCE_DIR}/applications/fastsum
      EXTRA_LINK fastsum_threads_lib fastsum_kernels)
  else()
    nfft_add_julia_module(fastsum libfastsumjulia.c
      EXTRA_INC ${PROJECT_SOURCE_DIR}/applications/fastsum
      EXTRA_LINK fastsum_lib fastsum_kernels)
  endif()
endif()
