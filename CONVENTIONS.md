# Code conventions used internally by NFFT3 (not in API)

## Common names

* `R`          : real type, (typically fftw_real)
* `E`          : real type for local variables (possibly extra precision)
* `C`          : complex type
* `A`          : assert
* `CK`         : check
* `X(...)`     : used for mangling of external names (see below)
* `Y(...)`     : used for mangling of external names (see below)
* `FFTW(...)`  : used for mangling of FFTW API names

## Name Mangling

Use Y(foo) for external names instead of nfft_foo.

Y(foo) expands to nfftf_foo, nfftl_foo, or nfft_foo, depending
on the precision.

X(foo) also expands to the corresponding prefix, but is local 
to each module. I.e. X(foo) will expand to nfct_foo in the 
NFCT module, but Y(foo) will always be nfft_foo.

FFTW(foo) expands to fftw_foo.

Names that are not exported do not need to be mangled.

## Indentation and Spacing

Indentation is processed using two spaces:
```
// Level 1.
  // Level 2.
```

Operators and operands, as well, should be padded with a single space
character such that every token is separated from neighbouring ones by a
space:
```
int a = 0x1 + 0x2;
```
   
