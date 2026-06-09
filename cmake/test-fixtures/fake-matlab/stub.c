/* Stub MATLAB C API — empty bodies so a mex links at test time (no runtime). */
#include <stddef.h>
double *mxGetPr(const void*p){return 0;} double *mxGetPi(const void*p){return 0;}
void *mxGetData(const void*p){return 0;} void *mxGetImagData(const void*p){return 0;}
size_t mxGetM(const void*p){return 0;} size_t mxGetN(const void*p){return 0;}
size_t mxGetNumberOfElements(const void*p){return 0;} size_t mxGetNumberOfDimensions(const void*p){return 0;}
double mxGetScalar(const void*p){return 0;} void *mxGetCell(const void*p,size_t i){return 0;}
int mxIsDouble(const void*p){return 0;} int mxIsSingle(const void*p){return 0;} int mxIsComplex(const void*p){return 0;}
int mxIsChar(const void*p){return 0;} int mxIsCell(const void*p){return 0;} int mxIsScalar(const void*p){return 0;}
int mxGetString(const void*p,char*s,size_t n){return 0;}
void *mxCreateDoubleMatrix(size_t m,size_t n,int c){return 0;}
void *mxCreateDoubleScalar(double d){return 0;}
void *mxCreateNumericArray(size_t n,const size_t*d,int cl,int cx){return 0;}
void *mxMalloc(size_t n){return 0;} void mxFree(void*p){}
void mexErrMsgTxt(const char*s){} void mexWarnMsgIdAndTxt(const char*a,const char*b,...){}
int mexPrintf(const char*f,...){return 0;} int mexEvalString(const char*s){return 0;}
void mexMakeMemoryPersistent(void*p){} int mexAtExit(void(*f)(void)){return 0;}
