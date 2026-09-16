import ctypes
import io
import os
import platform
import sys

from cpuinfo import get_cpu_info
from packaging.version import Version


# Create dummy classes for plans
class nfft_plan(ctypes.Structure):
    pass


class nfst_plan(ctypes.Structure):
    pass


class nfct_plan(ctypes.Structure):
    pass


class nfsft_plan(ctypes.Structure):
    pass


class fsft_plan(ctypes.Structure):
    pass


class fastsum_plan(ctypes.Structure):
    pass


glibcver = ""
is_macos = os.name != "nt" and os.uname().sysname == "Darwin"
# Apple Silicon Macs ship a single NEON build; they have no AVX/AVX2/SSE2
# instruction sets, so the x86 SIMD-level dispatch below does not apply.
is_arm64_macos = is_macos and platform.machine() in ("arm64", "aarch64")

# Determine the file extension for shared libraries based on the operating system
if os.name == "nt":  # Windows
    ending = ".dll"
elif is_macos:  # macOS
    ending = ".dylib"
else:  # Linux
    ending = ".so"
    glibcver = "glibc2.40"
    if Version(os.confstr("CS_GNU_LIBC_VERSION").split(" ")[1]) <= Version("2.35"):
        glibcver = "glibc2.22"

if is_arm64_macos:
    flag = "arm64"
elif "avx2" in get_cpu_info()["flags"]:
    flag = "AVX2"
elif "avx" in get_cpu_info()["flags"]:
    flag = "AVX"
elif "sse2" in get_cpu_info()["flags"]:
    flag = "SSE2"
else:
    raise RuntimeError("CPU type not supported")

# Check for CPU features and adjust library paths
package_dir = os.path.dirname(__file__)
lib_path_nfft = os.path.join(
    package_dir, "lib", flag, glibcver, "libnfftjulia" + ending
)
lib_path_nfct = os.path.join(
    package_dir, "lib", flag, glibcver, "libnfctjulia" + ending
)
lib_path_nfst = os.path.join(
    package_dir, "lib", flag, glibcver, "libnfstjulia" + ending
)
lib_path_nfsft = os.path.join(
    package_dir, "lib", flag, glibcver, "libnfsftjulia" + ending
)
lib_path_fastsum = os.path.join(
    package_dir, "lib", flag, glibcver, "libfastsumjulia" + ending
)

# Load the libraries
_nfftlib = ctypes.CDLL(lib_path_nfft)
_nfctlib = ctypes.CDLL(lib_path_nfct)
_nfstlib = ctypes.CDLL(lib_path_nfst)
_nfsftlib = ctypes.CDLL(lib_path_nfsft)
_fastsumlib = ctypes.CDLL(lib_path_fastsum)

# from .NFSFT import *
# from .FSFT import *
from .fastsum import *
from .flags import *
from .NFCT import *
# Import modules
from .NFFT import *
from .NFST import *

# Export functions and flags
__all__ = [
    "NFFT",
    "NFCT",
    "NFST",
    "NFSFT",
    "FSFT",
    "FASTSUM",
    "nfft_finalize_plan",
    "nfft_init",
    "nfft_trafo",
    "nfft_adjoint",
    "nfft_trafo_direct",
    "nfft_adjoint_direct",
    "nfct_finalize_plan",
    "nfct_init",
    "nfct_trafo",
    "nfct_adjoint",
    "nfct_transposed",
    "nfct_trafo_direct",
    "nfct_adjoint_direct",
    "nfct_transposed_direct",
    "nfst_finalize_plan",
    "nfst_init",
    "nfst_trafo",
    "nfst_adjoint",
    "nfst_trafo_direct",
    "nfst_adjoint_direct",
    "nfsft_finalize_plan",
    "nfsft_init",
    "nfsft_index",
    "nfsft_trafo",
    "nfsft_adjoint",
    "nfsft_trafo_direct",
    "nfsft_adjoint_direct",
    "nfsft_finalize_plan",
    "fsft_init",
    "fsft_index",
    "fsft_trafo",
    "fsft_adjoint",
    "fsft_trafo_direct",
    "fsft_adjoint_direct",
    "fastsum_finalize_plan",
    "fastsum_init",
    "fastsum_trafo",
    "fastsum_trafo_exact",
    "finalize_plan",
    "nffct_finalize_plan",
    "nffct_init",
    "nffct_trafo",
    "nffct_adjoint",
    "init",
    "trafo",
    "adjoint",
    "trafo_direct",
    "adjoint_direct",
    "trafo_exact",
    "PRE_PHI_HUT",
    "FG_PSI",
    "PRE_LIN_PSI",
    "PRE_FG_PSI",
    "PRE_PSI",
    "PRE_FULL_PSI",
    "MALLOC_X",
    "MALLOC_F_HAT",
    "MALLOC_F",
    "FFT_OUT_OF_PLACE",
    "FFTW_INIT",
    "PRE_ONE_PSI",
    "NFFT_SORT_NODES",
    "NFFT_OMP_BLOCKWISE_ADJOINT",
    "NFCT_SORT_NODES",
    "NFCT_OMP_BLOCKWISE_ADJOINT",
    "NFST_SORT_NODES",
    "NFST_OMP_BLOCKWISE_ADJOINT",
    "FFTW_MEASURE",
    "FFTW_DESTROY_INPUT",
    "FFTW_UNALIGNED",
    "FFTW_CONSERVE_MEMORY",
    "FFTW_EXHAUSTIVE",
    "FFTW_PRESERVE_INPUT",
    "FFTW_PATIENT",
    "FFTW_ESTIMATE",
    "FFTW_WISDOM_ONLY",
    "NFSFT_NORMALIZED",
    "NFSFT_USE_NDFT",
    "NFSFT_USE_DPT",
    "NFSFT_MALLOC_X",
    "NFSFT_MALLOC_F",
    "NFSFT_MALLOC_F_HAT",
    "NFSFT_PRESERVE_F_HAT",
    "NFSFT_PRESERVE_X",
    "NFSFT_PRESERVE_F",
    "NFSFT_DESTROY_F_HAT",
    "NFSFT_DESTROY_X",
    "NFSFT_DESTROY_F",
    "NFSFT_NO_DIRECT_ALGORITHM",
    "NFSFT_NO_FAST_ALGORITHM",
    "NFSFT_ZERO_F_HAT",
    "NFSFT_EQUISPACED",
    "f1_default_1d",
    "f1_default",
    "f2_default",
    "nfsft_default",
    "nfsft_nfft_default",
    "default_window_cut_off",
    "nfsft_default_nfft_cut_off",
    "nfsft_default_threshold",
]