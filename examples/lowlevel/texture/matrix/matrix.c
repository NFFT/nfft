#include<stdio.h>
#include<complex.h>
#include<stdlib.h>

#include<nfft3.h>
#include<texture_util.h>

int main()
{
	const char *grid_h_file = "h.in";
	const char *grid_r_file = "r.in";

	FILE *hf = fopen(grid_h_file, "r");
	FILE *rf = fopen(grid_r_file, "r");

	{
		int adjoint;
		int N, N1, N2;
		complex *omega, *x;
		double *h_phi, *h_theta, *r;
		texture_plan plan;

		read_grid(&N1, &N2, &h_phi, &h_theta, &r, hf, rf, stderr);

		scanf("%d%d", &adjoint, &N);

		omega = (complex *) smart_calloc(texture_flat_length(N), sizeof(complex));
		x = (complex *) smart_calloc(N1 * N2, sizeof(complex));

		texture_precompute(N);
		
		texture_init(&plan, N, N1, N2, omega, x, h_phi, h_theta, r);

		{
			int r;

			if (adjoint)
			{
				printf("%d %d\n", texture_flat_length(N), N1 * N2);

				for (r = 0; r < texture_flat_length(N); r++) {
					int s;

					omega[r] = 1;
					texture_set_omega(&plan, omega);

					texture_trafo(&plan);

					for(s = 0; s < N1 * N2; s++) {
						printf("%lg + %lgi ", creal(x[s]), -cimag(x[s]));
					}
					printf("\n");
					
					omega[r] = 0;
				}
			} else {
				printf("%d %d\n", N1 * N2, texture_flat_length(N));

				for (r = 0; r < N1 * N2; r++) {
					int s;
					
					x[r] = 1;
					texture_set_x(&plan, x);

					texture_adjoint(&plan);
					
					for (s = 0; s < texture_flat_length(N); s++) {
						printf("%lg + %lgi ", creal(omega[s]), -cimag(omega[s]));
					}
					printf("\n");

					x[r] = 0;
				}
			}
			
		}

		texture_finalize(&plan);

		texture_forget();

		free(omega);
		free(x);

		free(h_phi);
		free(h_theta);
		free(r);
	}

	fclose(hf);
	fclose(rf);

	return 0;
}
