#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <complex.h>
#include "nfft3.h"
#include "infft.h"

nfsft_plan* jnfsft_alloc(void) {
	nfsft_plan* p = nfft_malloc(sizeof(nfsft_plan));
	return p;
}

void jnfsft_init(nfsft_plan* p, int N, int M, unsigned int flags, unsigned int nfft_flags, int nfft_cutoff){
	nfsft_precompute(N,1000.0,0U,0U); // \todo: Make these variable for user
    nfsft_init_guru(p, N, M, flags, nfft_flags, nfft_cutoff);
}

double* jnfsft_set_x(nfsft_plan* p, double* X){
	int M = p->M_total;
	int j,c;
	nfsft_precompute_x(p);
	return p->x;
}

double _Complex* jnfsft_set_fhat(nfsft_plan* p, double _Complex* f_hat){
	int n = p->N_total;
	int k;
	for (k=0;k<n;k++)
		p->f_hat[k] = f_hat[k];
	return p->f_hat;
}

double _Complex* jnfsft_set_f(nfsft_plan* p, double _Complex* f){
	int M = p->M_total;
	int j;
	for (j=0;j<M;j++)
		p->f[j] = f[j];
	return p->f;
}

// nfsft trafo, return pointer to values for access by Julia if pointer isn't set
double _Complex* jnfsft_trafo(nfsft_plan* p){
	nfsft_trafo(p);
	return p->f;
}

// nfsft adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double _Complex* jnfsft_adjoint(nfsft_plan* p){
	nfsft_adjoint(p);
	return p->f_hat;
}

// nfsft trafo, return pointer to values for access by Julia if pointer isn't set
double _Complex* jnfsft_trafo_direct(nfsft_plan* p){
	nfsft_trafo_direct(p);
	return p->f;
}

// nfsft adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double _Complex* jnfsft_adjoint_direct(nfsft_plan* p){
	nfsft_adjoint_direct(p);
	return p->f_hat;
}

// nfsft plan finalizer
void jnfsft_finalize(nfsft_plan* p){
	nfsft_finalize(p);
	nfft_free(p);
	nfsft_forget();
}