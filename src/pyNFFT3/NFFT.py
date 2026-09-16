import ctypes

import numpy as np

from . import _nfftlib, nfft_plan
from .flags import *

# Set arugment and return types for functions
_nfftlib.jnfft_init.argtypes = [
    ctypes.POINTER(nfft_plan),
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.c_int32,
    ctypes.c_uint32,
    ctypes.c_uint32,
]

_nfftlib.jnfft_alloc.restype = ctypes.POINTER(nfft_plan)
_nfftlib.jnfft_finalize.argtypes = (ctypes.POINTER(nfft_plan),)

_nfftlib.jnfft_set_x.argtypes = [
    ctypes.POINTER(nfft_plan),
    np.ctypeslib.ndpointer(np.float64, flags="C"),
]
_nfftlib.jnfft_set_x.restype = ctypes.POINTER(ctypes.c_double)
_nfftlib.jnfft_set_f.argtypes = [
    ctypes.POINTER(nfft_plan),
    np.ctypeslib.ndpointer(np.complex128, ndim=1, flags="C"),
]
_nfftlib.jnfft_set_f.restype = ctypes.POINTER(ctypes.c_double)
_nfftlib.jnfft_set_fhat.argtypes = [
    ctypes.POINTER(nfft_plan),
    np.ctypeslib.ndpointer(np.complex128, ndim=1, flags="C"),
]
_nfftlib.jnfft_set_fhat.restype = ctypes.POINTER(ctypes.c_double)

_nfftlib.jnfft_trafo.argtypes = [ctypes.POINTER(nfft_plan)]
_nfftlib.jnfft_trafo.restype = ctypes.POINTER(ctypes.c_double)
_nfftlib.jnfft_adjoint.argtypes = [ctypes.POINTER(nfft_plan)]
_nfftlib.jnfft_adjoint.restype = ctypes.POINTER(ctypes.c_double)
_nfftlib.jnfft_trafo_direct.argtypes = [ctypes.POINTER(nfft_plan)]
_nfftlib.jnfft_trafo_direct.restype = ctypes.POINTER(ctypes.c_double)
_nfftlib.jnfft_adjoint_direct.argtypes = [ctypes.POINTER(nfft_plan)]
_nfftlib.jnfft_adjoint_direct.restype = ctypes.POINTER(ctypes.c_double)


class NFFT:
    """
    Class to perform non-equispaced fast Fourier transforms (NFFT)
    considering a **D**-dimensional trigonometric polynomial.
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
        self.M = M  # number of nodes
        self.n = n  # oversampling per dimension
        self.m = m  # window size
        self.D = len(N)  # dimensions

        if any(x <= 0 for x in N):
            raise ValueError(f"Invalid N: {N}. Argument must be a positive integer")

        if sum(x % 2 for x in N) != 0:
            raise ValueError(f"Invalid N: {N}. Argument must be an even integer")

        if M <= 0:
            raise ValueError(f"Invalid M: {M}. Argument must be a positive integer")

        if n is None:
            self.n = (2 ** (np.ceil(np.log(self.N) / np.log(2)) + 1)).astype("int32")

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

    def nfft_finalize_plan(self):
        """
        Finalizes an NFFT plan.
        This function does not have to be called by the user.
        """

        if not self.init_done:
            raise ValueError("NFFT not initialized.")

        if not self.finalized:
            _nfftlib.jnfft_finalize(self.plan)
            self.finalized = True

    def finalize_plan(self):
        """
        Alternate call for **nfft_finalize_plan()**
        """
        return self.nfft_finalize_plan()

    def nfft_init(self):
        """
        Initializes the NFFT plan in C.
        This function does not have to be called by the user.
        """
        # Convert N and n to numpy arrays for passing them to C
        Nv = np.array(self.N, dtype=np.int32)
        n = np.array(self.n, dtype=np.int32)

        # Call init for memory allocation
        ptr = _nfftlib.jnfft_alloc()

        # Set the pointer
        self.plan = ctypes.cast(ptr, ctypes.POINTER(nfft_plan))

        # Initialize values
        _nfftlib.jnfft_init(
            self.plan,
            ctypes.c_int32(self.D),
            ctypes.cast(Nv.ctypes.data, ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int32(self.M),
            ctypes.cast(n.ctypes.data, ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int32(self.m),
            self.f1,
            self.f2,
        )
        self.init_done = True

    def init(self):
        """
        Alternate call for **nfft_init()**
        """
        return self.nfft_init()

    @property
    def x(self) -> np.ndarray:
        return self._X

    @x.setter
    def x(self, value: np.ndarray):
        if value is not None:
            if not self.init_done:
                self.nfft_init()
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
                _nfftlib.jnfft_set_x(self.plan, value), shape=(self.M * self.D,)
            ).reshape(shape)

    @property
    def f(self) -> np.ndarray:
        return self._f

    @f.setter
    def f(self, value: np.ndarray):
        if value is not None:
            if not self.init_done:
                self.nfft_init()
            if self.finalized:
                raise RuntimeError("Plan already finalized")
            if not (
                isinstance(value, np.ndarray)
                and value.dtype == np.complex128
                and value.flags["C"]
            ):
                raise RuntimeError("f has to be C-continuous, numpy complex128 array")

            self._f = np.ctypeslib.as_array(
                _nfftlib.jnfft_set_f(self.plan, value), shape=(self.M * 2,)
            ).view(np.complex128)

    @property
    def fhat(self) -> np.ndarray:
        return self._fhat

    @fhat.setter
    def fhat(self, value: np.ndarray):
        if value is not None:

            if not self.init_done:
                self.nfft_init()

            if self.finalized:
                raise RuntimeError("Plan already finalized")

            if not (
                isinstance(value, np.ndarray)
                and value.dtype == np.complex128
                and value.flags["C"]
            ):
                raise RuntimeError(
                    "fhat has to be C-continuous, numpy complex128 array"
                )

            Ns = np.prod(self.N)

            self._fhat = np.ctypeslib.as_array(
                _nfftlib.jnfft_set_fhat(self.plan, value), shape=(Ns * 2,)
            ).view(np.complex128)

    @property
    def num_threads(self) -> int:
        return _nfftlib.nfft_get_num_threads()

    def nfft_trafo(self):
        """
        Computes the NDFT using the fast NFFT algorithm for the provided nodes in **x** and coefficients in **fhat**.
        """
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFFT already finalized")

        if not hasattr(self, "_fhat"):
            raise ValueError("fhat has not been set.")

        if not hasattr(self, "_X"):
            raise ValueError("x has not been set.")
        self._f = (
            np.ctypeslib.as_array(_nfftlib.jnfft_trafo(self.plan), shape=(self.M * 2,))
            .view(np.complex128)
            .copy()
        )

    def trafo(self):
        """
        Alternative call for **nfft_trafo()**
        """
        return self.nfft_trafo()

    def nfft_trafo_direct(self):
        """
        Computes the NDFT via naive matrix-vector multiplication for the provided nodes in **x** and coefficients in **fhat**.
        """
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFFT already finalized")

        if self._fhat is None:
            raise ValueError("fhat has not been set.")

        if self._X is None:
            raise ValueError("x has not been set.")

        self._f = (
            np.ctypeslib.as_array(
                _nfftlib.jnfft_trafo_direct(self.plan), shape=(self.M * 2,)
            )
            .view(np.complex128)
            .copy()
        )

    def trafo_direct(self):
        """
        Alternative call for **nfft_trafo_direct()**
        """
        return self.nfft_trafo_direct()

    def nfft_adjoint(self):
        """
        Computes the adjoint NDFT using the fast adjoint NFFT algorithm for the provided nodes in **x** and coefficients in **f**.
        """
        Ns = np.prod(self.N)
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFFT already finalized")

        if not hasattr(self, "_f"):
            raise ValueError("f has not been set.")

        if not hasattr(self, "_X"):
            raise ValueError("x has not been set.")

        self._fhat = (
            np.ctypeslib.as_array(_nfftlib.jnfft_adjoint(self.plan), shape=(Ns * 2,))
            .view(np.complex128)
            .copy()
        )

    def adjoint(self):
        """
        Alternative call for **nfft_adjoint()**
        """
        return self.nfft_adjoint()

    def nfft_adjoint_direct(self):
        """
        Computes the adjoint NDFT using naive matrix-vector multiplication for the provided nodes in **x** and coefficients in **f**.
        """
        Ns = np.prod(self.N)
        # Prevent bad stuff from happening
        if self.finalized:
            raise RuntimeError("NFFT already finalized")

        if not hasattr(self, "_f"):
            raise ValueError("f has not been set.")

        if not hasattr(self, "_X"):
            raise ValueError("x has not been set.")

        self._fhat = (
            np.ctypeslib.as_array(
                _nfftlib.jnfft_adjoint_direct(self.plan), shape=(Ns * 2,)
            )
            .view(np.complex128)
            .copy()
        )

    def adjoint_direct(self):
        """
        Alternative call for **nfft_adjoint_direct()**
        """
        return self.nfft_adjoint_direct()
