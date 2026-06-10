#ifndef MEX_H
#define MEX_H
#include "matrix.h"
void mexErrMsgTxt(const char*); void mexWarnMsgIdAndTxt(const char*, const char*, ...);
int mexPrintf(const char*, ...); int mexEvalString(const char*);
void mexMakeMemoryPersistent(void*); int mexAtExit(void(*)(void));
void mexFunction(int, mxArray**, int, const mxArray**);
#endif
