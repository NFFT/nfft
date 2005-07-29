#include <texture.h>
#include <util.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void print_error(complex* vec1, complex* vec2, int n) {
	double infNorm = 0, twoNorm = 0;
	int i;

	for(i = 0; i < n; i++) {
		double absDiff = cabs(vec1[i] - vec2[i]); 
		if(absDiff > infNorm) {
			infNorm = absDiff;
		}

		twoNorm += absDiff * absDiff;
	}
	twoNorm = sqrt(twoNorm);

	printf("error in twoNorm: %g in infNorm: %g\n", twoNorm, infNorm);
}

void simple_transform_test() {
	texture_plan my_plan;
	int N = 2;
	int N1 = (flatIndex(N, N, N) + 1)/2 + 3;
	int N2 = 1;
	int i;
	complex *f_hat_ref; 

	printf("simple transform test:\n");
	
	texture_init(&my_plan, N, N1, N2);

	f_hat_ref = (complex*) malloc(sizeof(complex) * my_plan.N_total);
	for(i = 0; i < my_plan.N_total; i++) {
		f_hat_ref[i] = i + 1;
		my_plan.f_hat[i] = f_hat_ref[i];
	}

	texture_trafo(&my_plan);

	vpr_complex(my_plan.f_hat, my_plan.N_total, "f_hat");
	vpr_complex(my_plan.f, my_plan.M_total, "f - the transform of f_hat");

	texture_adjoint(&my_plan);

	vpr_complex(my_plan.f_hat, my_plan.N_total, 
			"f_hat - the adjoint transform of f");

	free(f_hat_ref);
	texture_finalize(&my_plan);
}

void simple_solver_test() {
	texture_plan my_plan;
	itexture_plan my_iplan;
	int N = 2;
	int N1 = (flatIndex(N, N, N) + 1)/2 + 3;
	int N2 = 1;
	int i;
	complex *f_hat_ref; 

	printf("simple solver test:\n");
	
	texture_init(&my_plan, N, N1, N2);
	itexture_init(&my_iplan, &my_plan);

	f_hat_ref = (complex*) malloc(sizeof(complex) * my_plan.N_total);
	for(i = 0; i < my_plan.N_total; i++) {
		f_hat_ref[i] = i + 1;
		my_plan.f_hat[i] = f_hat_ref[i];
	}

	texture_trafo(&my_plan);
	vpr_complex(my_plan.f, my_plan.M_total, "f - the samples");
	vpr_complex(my_plan.f_hat, my_plan.N_total, 
			"f_hat - the corresponding frequencies");

	// destroying f_hat
	texture_adjoint(&my_plan);

	cp_complex(my_iplan.y, my_plan.f, my_plan.M_total);
	cp_complex(my_iplan.f_hat_iter, my_plan.f_hat, my_plan.N_total);
/*	for(i = 0; i < my_plan.N_total; i++) {
		my_iplan.f_hat_iter[i] = 0;
	}*/

	itexture_before_loop(&my_iplan);
	for(i = 0; i < 4; i++)	{
		vpr_complex(my_iplan.f_hat_iter, my_plan.N_total, "f_hat_iter - the guess");

		SWAP_complex(my_iplan.f_hat_iter, my_plan.f_hat);
		texture_trafo(&my_plan);
		print_error(my_iplan.y, my_plan.f, my_plan.M_total);
		SWAP_complex(my_iplan.f_hat_iter, my_plan.f_hat);

		itexture_loop_one_step(&my_iplan);
	}

	free(f_hat_ref);
	itexture_finalize(&my_iplan);
	texture_finalize(&my_plan);
}

int main() {
	simple_transform_test();
	simple_solver_test();
	return 0;
}
