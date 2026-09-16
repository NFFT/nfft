import ctypes

import numpy as np

from . import _nfstlib, nfst_plan
from .flags import *

# Set arugment and return types for functions
_nfstlib.jnfst_init.argtypes = [
    ctypes.POINTER(nfst_plan),
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_int32,
    ctypes.c_uint32,
    ctypes.c_uint32,
]

_nfstlib.jnfst_alloc.restype = ctypes.POINTER(nfst_plan)
_nfstlib.jnfst_finalize.argtypes = (ctypes.POINTER(nfst_plan),)

_nfstlib.jnfst_set_x.argtypes = [
    ctypes.POINTER(nfst_plan),
    np.ctypeslib.ndpointer(np.float64, flags="C"),
]
_nfstlib.jnfst_set_x.restype = ctypes.POINTER(ctypes.c_double)
_nfstlib.jnfst_set_f.argtypes = [
    ctypes.POINTER(nfst_plan),
    np.ctypeslib.ndpointer(np.float64, ndim=1, flags="C"),
]
_nfstlib.jnfst_set_f.restype = ctypes.POINTER(ctypes.c_double)
_nfstlib.jnfst_set_fhat.argtypes = [
    ctypes.POINTER(nfst_plan),
    np.ctypeslib.ndpointer(np.float64, ndim=1, flags="C"),
]
_nfstlib.jnfst_set_fhat.restype = ctypes.POINTER(ctypes.c_double)

_nfstlib.jnfst_trafo.argtypes = [ctypes.POINTER(nfst_plan)]
_nfstlib.jnfst_trafo.restype = ctypes.POINTER(ctypes.c_double)
_nfstlib.jnfst_adjoint.argtypes = [ctypes.POINTER(nfst_plan)]
_nfstlib.jnfst_adjoint.restype = ctypes.POINTER(ctypes.c_double)
_nfstlib.jnfst_trafo_direct.argtypes = [ctypes.POINTER(nfst_plan)]
_nfstlib.jnfst_trafo_direct.restype = ctypes.POINTER(ctypes.c_double)
_nfstlib.jnfst_adjoint_direct.argtypes = [ctypes.POINTER(nfst_plan)]
_nfstlib.jnfst_adjoint_direct.restype = ctypes.POINTER(ctypes.c_double)


class NFST:
    """
    Class to perform non-equispaced fast sine transforms (NFST)
    With a dimension of **D**.
    Just **N** and **M** are required for initializing a plan.
    """

    def __init__(
        self,
        N: np.ndarray,
        M: int,
        n: np.ndarray = None,
        m: int = default_window_cut_off,
        f1: ctypes.c_uint32 = None,
        f2: ctypes.c_uint32 = f2_default,
    ):
        self.plan = None
        self.N = N  # bandwidth tuple
        N_ct = N.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        self.M = M  # number of nodes
        self.n = n  # oversampling per dimension
        self.m = m  # window size
        self.D = len(N)  # dimensions

        if any(x <= 0 for x in N):
            raise ValueError(f"Invalid N: {N}. Argument must be a positive integer")

        if M <= 0:
            raise ValueError(f"Invalid M: {M}. Argument must be a positive integer")

        if n is None:
            self.n = (2 ** (np.ceil(np.log(self.N) / np.log(2)) + 1)).astype("int32")
            n_ct = self.n.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))

        if any(x <= 0 for x in self.n):
            raise ValueError(
                f"Invalid n: {self.n}. Argument must be a positive integer"
            )

        if any(x <= y for x, y in zip(self.n, N)):
            raise ValueError(f"Invalid n: {self.n}. Argument must fulfil n_i > N_i")

        if sum(x % 2 for x in self.n) != 0:
            raise ValueError(f"Invalid n: {self.n}. Argument must be an even integer")

        if m <= 0:
            raise ValueError(f"Invalid m: {m}. Argument must be a positive integer")

        if f1 is None:
            self.f1 = f1_default if self.D > 1 else f1_default_1d
        else:
            self.f1 = f1

        self.f2 = f2  # FFTW flags
        self.init_done = False  # bool for plan init
        self.finalized = False  # bool for finalizer
        self._X = None  # nodes, will be set later
        self._f = None  # function values
        self._fhat = None  # Fourier coefficients

    def __del__(self):
        self.finalize_plan()

    def nfst_finalize_plan(self):
        """
        Finalizes an NFST plan.
        This function does not have to be called by the user.
        """

        if not self.init_done:
            raise ValueError("NFST not initialized.")

        if not self.finalized:
            _nfstlib.jnfst_finalize(self.plan)
            self.finalized = True

    def finalize_plan(self):
        """
        Alternative call for **nfst_finalize_plan()**
        """
        return self.nfst_finalize_plan()

    def nfst_init(self):
        """
        Initializes the NFST plan in C.
        This function does not have to be called by the user.
        """
        # Convert N and n to numpy arrays for passing them to C
        Nv = np.array(self.N, dtype=np.int32)
        n = np.array(self.n, dtype=np.int32)

        # Call init for memory allocation
        ptr = _nfstlib.jnfst_alloc()

        # Set the pointer
        self.plan = ctypes.cast(ptr, ctypes.POINTER(nfst_plan))

        # Initialize values
        _nfstlib.jnfst_init(
            self.plan,
            ctypes.c_int(self.D),
            ctypes.cast(Nv.ctypes.data, ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int(self.M),
            ctypes.cast(n.ctypes.data, ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int(self.m),
            self.f1,
            self.f2,
        )
        self.init_done = True

    def init(self):
        """
        Alternative call for **nfst_init()**
        """
        return self.nfst_init()

    @property
    def x(self) -> np.ndarray:
        return self._X

    @x.setter
    def x(self, value: np.ndarray):
        if value is not None:
            if not self.init_done:
                self.nfst_init()
            if self.finalized:
                raise RuntimeError("Plan already finalized")
            if not (
                isinstance(value, np.ndarray)
                and value.dtype == np.float64
                and value.flags["C"]
            ):
                raise RuntimeError("x has to be C-continuous, numpy float64 array")
            if self.D == 1:
                shape = self.M
            else:
                shape = (self.M, self.D)
            self._X = np.ctypeslib.as_array(
                _nfstlib.jnfst_set_x(self.plan, value), shape=(self.M * self.D,)
            ).reshape(shape)

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self, value: np.ndarray):
        if value is not None:
            if not self.init_done:
                self.nfst_init()
            if self.finalized:
                raise RuntimeError("Plan already finalized")
            if not (
                isinstance(value, np.ndarray)
                and value.dtype == np.float64
                and value.flags["C"]
            ):
                raise RuntimeError("f has to be C-continuous, numpy float64 array")
            self._f = np.ctypeslib.as_array(
                _nfstlib.jnfst_set_f(self.plan, value), shape=(self.M,)
            )

    @property
    def fhat(self):
        return self._fhat

    @fhat.setter
    def fhat(self, value: np.ndarray):
        if value is not None:
            Ns = np.prod(self.N - 1)
            if not self.init_done:
                self.nfst_init()
            if self.finalized:
                raise RuntimeError("Plan already finalized")
            if not (isinstance(value, np.ndarray) and value.dtype == np.float64):
                raise RuntimeError("fhat has to be a numpy float64 array")
            if not value.flags["C"]:
                raise RuntimeError("fhat has to be C-continuous")
            if value.size != Ns:
                raise ValueError(f"fhat has to be an array of size {Ns}")
            self._fhat = np.ctypeslib.as_array(
                _nfstlib.jnfst_set_fhat(self.plan, value), shape=(Ns,)
            )

    @property
    def num_threads(self) -> int:
        return _nfstlib.nfft_get_num_threads()

    def nfst_trafo(self):
        """
        Computes the NDFT via the fast NFST algorithm for the provided nodes in **x** and coefficients in **fhat**.
        """
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFST already finalized")

        if not hasattr(self, "_fhat"):
            raise ValueError("fhat has not been set.")

        if not hasattr(self, "_X"):
            raise ValueError("x has not been set.")
        self._f = np.ctypeslib.as_array(
            _nfstlib.jnfst_trafo(self.plan), shape=(self.M,)
        ).copy()

    def trafo(self):
        """
        Alternative call for **nfst_trafo()**
        """
        return self.nfst_trafo()

    def nfst_trafo_direct(self):
        """
        Computes the NDST via naive matrix-vector multiplication for provided nodes in **x** and coefficients in **fhat**.
        """
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFST already finalized")

        if self._fhat is None:
            raise ValueError("fhat has not been set.")

        if self._X is None:
            raise ValueError("x has not been set.")
        self._f = np.ctypeslib.as_array(
            _nfstlib.jnfst_trafo_direct(self.plan), shape=(self.M,)
        ).copy()

    def trafo_direct(self):
        """
        Alternative call for nfst_trafo_direct()
        """
        return self.nfst_trafo_direct()

    def nfst_transposed(self):
        """
        Computes the transposed NDST via the fast transposed NFST algorithm for the provided nodes in **x** and coefficients in **f**.
        """
        Ns = np.prod(self.N - 1)
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFST already finalized")

        if self._f is None:
            raise ValueError("f has not been set.")

        if self._X is None:
            raise ValueError("x has not been set.")
        self._fhat = np.ctypeslib.as_array(
            _nfstlib.jnfst_adjoint(self.plan), shape=(Ns,)
        ).copy()

    def nfst_transposed_direct(self):
        """
        Computes the transposed NDST via naive matrix-vector multiplication for provided nodes in **x** and coefficients in **f**.
        """
        Ns = np.prod(self.N - 1)
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFST already finalized")

        if self._f is None:
            raise ValueError("f has not been set.")

        if self._X is None:
            raise ValueError("x has not been set.")
        self._fhat = np.ctypeslib.as_array(
            _nfstlib.jnfst_adjoint_direct(self.plan), shape=(Ns,)
        ).copy()

    def nfst_adjoint(self):
        """
        Alternative call for nfst_transposed()
        """
        return self.nfst_transposed()

    def nfst_adjoint_direct(self):
        """
        Alternative call for nfst_transposed_direct()
        """
        return self.nfst_transposed_direct()

    def adjoint(self):
        """
        Alternative call for nfst_adjoint()
        """
        return self.nfst_adjoint()

    def adjoint_direct(self):
        """
        Alternative call for nfst_adjoint_direct()
        """
        return self.nfst_adjoint_direct()
