#ifndef MATRIX_H
#define MATRIX_H
#include <stddef.h>
typedef struct mxArray_tag mxArray;
typedef size_t mwSize; typedef size_t mwIndex;
typedef enum { mxUNKNOWN_CLASS=0, mxDOUBLE_CLASS, mxSINGLE_CLASS } mxClassID;
typedef enum { mxREAL=0, mxCOMPLEX } mxComplexity;
double *mxGetPr(const mxArray*); double *mxGetPi(const mxArray*);
void *mxGetData(const mxArray*); void *mxGetImagData(const mxArray*);
size_t mxGetM(const mxArray*); size_t mxGetN(const mxArray*);
size_t mxGetNumberOfElements(const mxArray*); size_t mxGetNumberOfDimensions(const mxArray*);
double mxGetScalar(const mxArray*); mxArray *mxGetCell(const mxArray*, size_t);
int mxIsDouble(const mxArray*); int mxIsSingle(const mxArray*); int mxIsComplex(const mxArray*);
int mxIsChar(const mxArray*); int mxIsCell(const mxArray*); int mxIsScalar(const mxArray*);
int mxGetString(const mxArray*, char*, size_t);
mxArray *mxCreateDoubleMatrix(size_t, size_t, mxComplexity);
mxArray *mxCreateDoubleScalar(double);
mxArray *mxCreateNumericArray(size_t, const size_t*, mxClassID, mxComplexity);
void *mxMalloc(size_t); void mxFree(void*);
#endif
