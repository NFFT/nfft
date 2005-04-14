#include "nfft.h"
#include "utils.h"

void short_nfft_trafo(nfft_plan* this_plan, nfft_plan* plan1d);
void short_nfft_trafo_horner(nfft_plan* this_plan, nfft_plan* plan_1d, fftw_complex* exp_omega1);

void short_nfft_trafo_3d_1(nfft_plan* this_plan, nfft_plan* plan1d);
void short_nfft_trafo_3d_2(nfft_plan* this_plan, nfft_plan* plan2d);
