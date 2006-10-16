#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>

#include<texture_util.h>

#include"data_util.h"

int N1, N2;
double *h_phi, *h_theta, *r;
complex *x;

int
main ()
{
  FILE *fid = fopen ("grid_size", "r");
  int x1, x2, x3;
  int i, j;

  fscanf (fid, "%d%d", &N1, &N2);
  fclose (fid);

  h_phi = smart_malloc (N1 * sizeof (double));
  h_theta = smart_malloc (N1 * sizeof (double));
  r = smart_malloc (2 * N1 * N2 * sizeof (double));
  x = smart_malloc (N1 * N2 * sizeof (complex));

  fid = fopen ("h_files", "r");
  for (i = 0; i < N1; i++)
    {
      double norm;
      fscanf (fid, "%*[^0-9]%1d%1d%1d", &x1, &x2, &x3);
      norm = sqrt (x1 * x1 + x2 * x2 + x3 * x3);

      h_theta[i] = normalise_theta (acos (x3 / norm));
      h_phi[i] = normalise_phi (atan2 (x2, x1));
    }
  fclose (fid);

  fid = fopen ("h_files", "r");
  for (i = 0; i < N1; i++)
    {
      FILE *pfid;
      char pfname[20];
      char pfpath[40];

      fscanf (fid, "%s\n", pfname);
      sprintf (pfpath, "data/%s", pfname);
      pfid = fopen (pfpath, "r");

      for (j = 0; j < N2; j++)
	{
	  double dummy;
	  fscanf (pfid, "%lg%lg%lg", &(r[2 * (i * N2 + j) + 1]),
		  &(r[2 * (i * N2 + j)]), &dummy);
	  x[i * N2 + j] = dummy;
	}
      fclose (pfid);
    }
  fclose (fid);
  normalise_r (N2, r);

  fid = fopen ("h_file", "w");
  fprintf (fid, "Polefigures\n");
  fprintf (fid, "# Data set for the publication.\n");
  write_h (N1, h_phi, h_theta, fid);
  fclose (fid);

  // The result of the following piece of code is that all nodesets are equal.
  /*
     for (i = 0; i < N1; i++)
     {
     char name[20];

     sprintf (name, "r_file%d", i);

     fid = fopen (name, "w");
     write_r (N2, &(r[2 * i * N2]), fid);
     fclose (fid);
     }
   */

  fid = fopen ("r_file", "w");
  fprintf (fid, "Nodes\n");
  fprintf (fid, "# Data set for the publication.\n");
  write_r (N2, r, fid);
  fclose (fid);

  fid = fopen ("samples", "w");
  fprintf (fid, "Samples\n");
  fprintf (fid, "# Data set for the publication.\n");
  write_x (N1, N2, x, fid);
  fclose (fid);

  free (h_phi);
  free (h_theta);
  free (r);
  free (x);

  return 0;
}
