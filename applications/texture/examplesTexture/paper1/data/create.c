#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include<util.h>
#include<texture_util.h>
#include<data_util.h>

int N1, N2;
int new_N1, new_N2;
double *h_phi, *h_theta, *r;
const char *input_file_prefix = "raw/";
const char *output_file_prefix = "";
char input_file_name[100], output_file_name[100];
char input_file_path[200], output_file_path[200];
char file_type;
int mirror;

void input()
{
	// input filename
	printf("Name of the input file: ");
	scanf("%s", input_file_name);

	// pole figures or nodes
	do {
		printf("Pole figures (h) or nodes (r): ");
		scanf(" %c", &file_type);
	} while (file_type != 'h' && file_type != 'r');

	// mirror?
	printf("Mirror data: ");
	scanf("%d", &mirror);

	// size 
	printf("Size (N1 / N2): ");
	scanf("%d", &N1);
	new_N1 = (mirror) ? 2 * N1 : N1;
	N2 = N1;
	new_N2 = new_N1;

	h_phi = smart_malloc(new_N1 * sizeof(double));
	h_theta = smart_malloc(new_N1 * sizeof(double));
	r = smart_malloc(new_N2 * 2 * sizeof(double));

	// output filename
	printf("Name of the output file: ");
	scanf("%s", output_file_name);
}

void convert()
{
	FILE *input;

	sprintf(input_file_path, "%s%s", input_file_prefix, input_file_name);
	input = fopen(input_file_path, "r");

	if (file_type == 'h') {
		int i;

		for (i = 0; i < N1; i++) {
			fscanf(input, "%lg %lg", &(h_theta[i]), &(h_phi[i]));
			if (PI < h_theta[i] && h_theta[i] <= 3.1415927) {
				h_theta[i] = PI;
			}
		}

		if (mirror) {
			expand_h(new_N1, h_phi, h_theta);
		}
		
		normalise_h(new_N1, h_phi, h_theta);

	} else {
		int i;

		for (i = 0; i < N2; i++) {
			fscanf(input, "%lg %lg", &(r[2 * i + 1]), &(r[2 * i]));
			if (PI < r[2 * i + 1] && r[2 * i + 1] <= 3.1415927) {
				r[2 * i + 1] = PI;
			}
		}

		if (mirror) {
			expand_r(new_N2, r);
		}

		normalise_r(new_N2, r);

	}

	fclose(input);
}

void output()
{
	FILE *output;

	sprintf(output_file_path, "%s%s", output_file_prefix, output_file_name);
	output = fopen(output_file_path, "w");

	if (file_type == 'h') {
		fprintf(output, "Polefigures\n");
		fprintf(output, "# Approximate equidistribution on the hemisphere.\n");
		fprintf(output, "# From file %s.\n", input_file_name);
		write_h(new_N1, h_phi, h_theta, output);
	} else {
		fprintf(output, "Nodes\n");
		fprintf(output, "# Approximate equidistribution on the sphere.\n");
		fprintf(output, "# From file %s.\n", input_file_name);
		write_r(new_N2, r, output);
	}

	fclose(output);
}

void cleanup()
{
	free(h_phi);
	free(h_theta);
	free(r);
}

int main()
{
	input();

	convert();

	output();

	cleanup();

	return 0;
}
