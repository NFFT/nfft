#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include "util.h"
#include "nfft3.h"



/*Struktur fuer die Knoten des Octrees*/
struct octreenode{
    struct octreenode *child[8]; /*Zeiger auf die Kinder*/
    double cell_center[3]; /* Mittelpunkte der Wuerfel*/
    double cell_vertices[8][3]; /*Eckpunkte der Wuerfel*/
    double atom_i[3];  /* Array um das Atom, das in diesem Wuerfel liegt einzufuegen*/
    int atom; /* wenn atom =1 wurde Atom eingefuegt*/
    double sidelength; /*Kantenlaenge des Wuerfels*/
    int sas; /*Zuordnung ob Knoten zum SAS-Volume(ses=1) zum SES-Surface(sas=2) oder außerhalb des Moleküls liegt (sas=0)*/
    double radius_i; /*Radius des eingefuegten Atoms*/
    };

/* Initialisierung von refine_tree(), traverse_tree, abstand und free_tree */
struct octreenode* refine_tree(struct octreenode *this,double c_a[], double l);
void traverse_tree(struct octreenode *this, double *atomkoordinaten, double atomradius[], int molekuelzahl);
double abstand(double vektor1[], double vektor2[]);
void free_tree(struct octreenode *this);

struct b_l
{
  double n;                                /* expansion degree */
  double p;                             /* period */
  complex decay;                             /* rate of decay parameter*/
  double R_i;                            /* van der Waals radius*/
  complex *b;                           /* expansion coefficients */
};

int main(int argc, char *argv[ ])
{
    int i;
    int j;
    int k,l,r;
    double s,t,u;

    /*Atomzahl von Molekuel A und B einlesen und umwandeln in int*/
    /*Anzahl der Atome der beiden Molekuele; diese entspricht der Anzahl der Atome des
    groesseren Molekuels*/
   int atomzahl_b = atoi(argv[2]);
   int atomzahl_a_orig = atoi(argv[1]);
   int atomzahl_a= atomzahl_a_orig+ 2500;

    /*Kantenlaenge des Wuerfels*/
    double kantenlaenge;


    /*Atomkoordinaten der beiden Molekuele*/
    double *molekuel_a;
   double *molekuel_a_extra;
    double *molekuel_b;
    double *molekuel_b_neu;

   /*Speicherplatz für die Atomkoordinaten fuer Molekuel A und B*/

    molekuel_a = (double *)malloc(3*atomzahl_a_orig*sizeof(double));
    molekuel_b = (double *)malloc(3*atomzahl_b*sizeof(double));
    molekuel_a_extra = (double *)malloc(3*3000*sizeof(double));


    /*Arrays fuer die Radien der Atome*/
    char atomrad_mom[atomzahl_a_orig];
    double atomradien[atomzahl_a_orig];

    /*Array fuer die SAS-Werte des Molekuels B*/
    double *sas_B;

    /*Array fuer die gamma-Werte von Molekuel A und B*/
    double *gamma_A;
    double *gamma_B;


    /*file mit den SAS-Werten von Molekuel B*/
    FILE *saswerte_B;
    /*file mit den Atomkoordinaten von Molekuel A und B*/
    FILE *atomkoord_A;
    FILE *atomkoord_B;
    FILE *atomkoord_A_extra;
    /*FILE *atomkoord_A_out;
    FILE *atomkoord_B_out;*/
    FILE *atomrad_A;

   /*Array fuer das Zentrum der Molekuele*/
   double zentrum_A[3];
   double zentrum_B[3];

    /*Arrays fuer die Summe der Atomkoordinaten*/
   double sum_A[3];
   double sum_B[3];

    /*parameter p (entspricht hier max_Abstand_A und max_Abstand_B)to embed particle data in volumetric interval I und Variablen fuer
   maximalen Abstand zwischen den Atomkoordinaten*/
   double abstand_mom;
   double max_abstand_A=0.0;
   double max_abstand_B=0.0;
   double p;

   /**/
   double vektor_1[3];
   double vektor_2[3];

   /*Initialisierung der Wurzel des octrees*/
   struct octreenode *root = NULL;

   /*van der Waals-Radius R_k einmal fest mit dem Wert 1 und dann jeweils fuer die in AS vorkommenden Atomen*/
   const double van_der_waals_radius=1.0;
   const double vdW_radius_hydrogen=1.2;
   const double vdW_radius_carbon=1.7;
   const double vdW_radius_oxygen=1.4;
   const double vdW_radius_nitrogen=1.55;
   const double vdW_radius_sulphur=1.85;
   const double vdw_radius_water= 1.8;

   /*rotation angles*/
   double alpha;
   double beta;
   double gamma;

    /*parameters for NFFT*/
    nfft_plan plan_a;
    nfft_plan plan_b;
    nfft_plan plan_a_b;
    int M;

    complex *b_l_3d;
    complex *h_k;

    struct b_l this;
    this.n = 64.0*64.0*64.0;
    this.p = 1.0;
    this.decay =-0.1+I*0.0;
    this.R_i = 1.0;
    this.b = (complex*)fftw_malloc(this.n*sizeof(complex));


    /*Oeffnen der Dateien*/
    saswerte_B = fopen(argv[6],"r");
    atomkoord_A = fopen(argv[3],"r");
    atomkoord_B = fopen(argv[4],"r");
    atomkoord_A_extra=fopen(argv[5],"r");
    atomrad_A=fopen(argv[7],"r");

     /* Einlesen der Atomkoordinaten fuer Molekuel A */
   if(atomkoord_A == NULL || atomrad_A==NULL)
   {
      printf("Dateien konnte nicht geoeffnet werden!\n");
   }
   else
   {
      for(i=0;i<atomzahl_a_orig;i++)
        {
           fscanf(atomkoord_A, "%lf %lf %lf \n",
	   &molekuel_a[3*i+0],&molekuel_a[3*i+1],&molekuel_a[3*i+2]);

	   fscanf(atomrad_A, "%c \n",&atomrad_mom[i]);
	   if(atomrad_mom[i]== 'H')
	   {
               atomradien[i] = vdW_radius_hydrogen;
	   }
	   if(atomrad_mom[i]== 'C')
	   {
	       atomradien[i] = vdW_radius_carbon;
	   }
	   if(atomrad_mom[i]== 'O')
	   {
	       atomradien[i] = vdW_radius_oxygen;
	   }
	   if(atomrad_mom[i]== 'N')
	   {
	       atomradien[i] = vdW_radius_nitrogen;
	   }
	   if(atomrad_mom[i]== 'S')
	   {
	       atomradien[i] = vdW_radius_sulphur;
	   }
        }
   }


   /* for(i=0;i<atomzahl_a_orig;i++)
   {
        for(j=0;j<3;j++)
        {
            printf("%lf ",molekuel_a[3*i+j]);
        }
        printf("\n");
   }


   for(i=0;i<atomzahl_a_orig;i++)
      {
          printf("%lf ", atomradien[i]);
          printf("\n");
      }
    */

    /*Schliessen der Datei*/
    fclose(atomrad_A);

    /*Zentrum des Molekuels A*/
     for(i=0;i<atomzahl_a_orig; i++)
      for( j=0; j<3; j++)
         {
	    sum_A[j] += molekuel_a[3*i+j];
	 }
     for( j=0; j<3; j++)
     {
        zentrum_A[j] = ((sum_A[j])/(atomzahl_a_orig));
     }

     /*Parameter P*/
     for(i=0;i<atomzahl_a_orig;i++)
         for(j=0;j<atomzahl_a_orig;j++)
         {
	     abstand_mom=sqrt(pow((molekuel_a[3*i+0]-molekuel_a[3*j+0]),2)
	                 +pow((molekuel_a[3*i+1]-molekuel_a[3*j+1]),2)
	                 +pow((molekuel_a[3*i+2]-molekuel_a[3*j+2]),2));

	    if(max_abstand_A < abstand_mom)
	       {
	         max_abstand_A = abstand_mom;
	       }
	 }
    p = max_abstand_A;


    /*memory allocation for root*/
    root = (struct octreenode *) malloc(sizeof (root));

    kantenlaenge = p + 4* vdw_radius_water +2; /*Kantenlaenge des Wuerfels ist maximaler Abstand + 4mal der Radius von H2O */

    /*Aufruf von refine_tree und traverse_tree*/
   /* refine_tree(root,zentrum_A,kantenlaenge);
    traverse_tree(root, molekuel_a,atomradien ,atomzahl_a_orig);
    free_tree(root);*/


    /* Einlesen der Atomkoordinaten fuer Molekuel A und B */
    if(atomkoord_A_extra == NULL)
    {
       printf("Datei A konnte nicht geoeffnet werden!\n");
    }
    else
    {
        /*for(i=0;i<atomzahl_a_orig;i++)
	{
	   fscanf(atomkoord_A, "%lf %lf %lf \n", &molekuel_a[3*i+0 ],&molekuel_a[3*i+1],&molekuel_a[3*i+2]);
	}
	for(i=atomzahl_a_orig;i<atomzahl_a;i++)
	{
             fscanf(atomkoord_A_extra, "%lf %lf %lf \n", &molekuel_a[3*i+0 ],&molekuel_a[3*i+1],&molekuel_a[3*i+2]);
	}*/

        i=0;
	 while(!feof(atomkoord_A_extra)&& i<=3000)
	 {

	   fscanf(atomkoord_A_extra, "%lf %lf %lf \n", &molekuel_a_extra[3*i+0 ],&molekuel_a_extra[3*i+1],&molekuel_a_extra[3*i+2]);

            if(i!=0)
	    {
	      abstand_mom = sqrt(pow((molekuel_a_extra[3*(i-1)+0]-molekuel_a_extra[3*i+0]),2)
	                       + pow((molekuel_a_extra[3*(i-1)+1]-molekuel_a_extra[3*i+1]),2)
	                       + pow((molekuel_a_extra[3*(i-1)+2]-molekuel_a_extra[3*i+2]),2));
	          if(abstand_mom<vdw_radius_water)
	          {
                      i--;
	          }
	       }
	    i++;
         }
	 atomzahl_a=(atomzahl_a_orig+i);
	 /*printf("%i ",atomzahl_a);*/
      }


      molekuel_a = (double *)realloc(molekuel_a,3*(atomzahl_a)*sizeof(double));
      molekuel_b_neu = (double *)malloc(3*atomzahl_a*sizeof(double));

      for(i=0;i<(atomzahl_a-atomzahl_a_orig);i++)
      {
         for(j=0;j<3;j++)
	 {
	    molekuel_a[3*(i+(atomzahl_a_orig))+j] = molekuel_a_extra[3*i+j];
	 }
      }



     /* for(i=0;i<atomzahl_a;i++)
      {
        for(j=0;j<3;j++)
        {
            printf("%lf ",molekuel_a[3*i+j]);
        }
        printf("\n");
      }
     */


    sas_B=(double *)malloc(3*atomzahl_b*sizeof(double));
    gamma_B=(double *)malloc(3*atomzahl_a*sizeof(double));
    gamma_A=(double *)malloc(3*atomzahl_a*sizeof(double));

    if(atomkoord_B == NULL)
    {
      printf("Datei B konnte nicht geoeffnet werden!\n");
    }
    else
    {
        for(i=0;i<atomzahl_b;i++)
        {
           fscanf(atomkoord_B, "%lf %lf %lf \n", &molekuel_b[3*i+0 ],&molekuel_b[3*i+1],&molekuel_b[3*i+2]);
        }
    }



     /* Einlesen der SAS-Werte fuer Molekuel B */
    if(saswerte_B == NULL)
    {
       printf("Datei mit SAS-Werten konnte nicht geoeffnet werden!\n");
    }
    else
    {
       for(i=0;i<atomzahl_b;i++)
         {
            fscanf(saswerte_B, "%lf \n",&sas_B[i]);
	    if(sas_B[i]>1)
	    {
               gamma_B[2*i+0]=1.0 ;
	       gamma_B[2*i+1]=0.0;
	    }
	    else
	    {
                gamma_B[2*i+0]=0.0 ;
	       gamma_B[2*i+1]=-9.0;
	    }
         }
	for(i=atomzahl_b;i<atomzahl_a;i++)
	{
            gamma_B[2*i+1]=-9.0;
	    gamma_B[2*i+0]=0.0;

	}
     }

     for(i=0;i<atomzahl_a_orig;i++)
     {
        gamma_A[2*i+1]=9.0;
	gamma_A[2*i+0]=0.0;
     }

     for(i=atomzahl_a_orig; i< atomzahl_a;i++)
     {
        gamma_A[2*i+0]=1.0;
	gamma_A[2*i+1]=0.0;
     }

      /*Schliessen der beiden Dateien*/
    fclose(atomkoord_A);
    fclose(atomkoord_B);

    for(j=0;j<3;j++)
    {
        sum_A[j]=0.0;
	sum_B[j]=0.0;
	zentrum_A[j]=0.0;
	zentrum_A[j]=0.0;
    }

   /*Zentrum der Molekuele A und B*/
   for(i=0;i<atomzahl_a; i++)
      for( j=0; j<3; j++)
         {
	    sum_A[j] += molekuel_a[3*i+j];
	 }

   for(i=0;i<atomzahl_b; i++)
      for( j=0; j<3; j++)
         {
	    sum_B[j] += molekuel_b[3*i+j];
	 }

   for( i=0; i<3; i++)
      {
         zentrum_A[i] = ((sum_A[i])/(atomzahl_a));
	 zentrum_B[i] = ((sum_B[i])/(atomzahl_b));
      }

   /*Parameter P*/
   for(i=0;i<atomzahl_a;i++)
      for(j=0;j<atomzahl_a;j++)
         {
	     abstand_mom=sqrt(pow((molekuel_a[3*i+0]-molekuel_a[3*j+0]),2)
	                 +pow((molekuel_a[3*i+1]-molekuel_a[3*j+1]),2)
	                 +pow((molekuel_a[3*i+2]-molekuel_a[3*j+2]),2));

	     if(max_abstand_A < abstand_mom)
	       {
	         max_abstand_A = abstand_mom;
	       }
	  }
     for(i=0;i<atomzahl_b;i++)
       for(j=0;j<atomzahl_b;j++)
          {
	     abstand_mom=sqrt(pow((molekuel_b[3*i+0]-molekuel_b[3*j+0]),2)
	                 +pow((molekuel_b[3*i+1]-molekuel_b[3*j+1]),2)
	                 +pow((molekuel_b[3*i+2]-molekuel_b[3*j+2]),2));

	     if(max_abstand_B < abstand_mom)
	        {
	           max_abstand_B = abstand_mom;
	        }
	  }

   /*Recenter und Rescale der particle center*/
   for(i=0;i<atomzahl_a;i++)
      for(j=0;j<3;j++)
         {
	   molekuel_a[3*i+j]= ((molekuel_a[3*i+j]-zentrum_A[j])/(2*max_abstand_A));

	 }

   for(i=0;i<atomzahl_b;i++)
      for(j=0;j<3;j++)
         {

	   molekuel_b[3*i+j]= ((molekuel_b[3*i+j]-zentrum_B[j])/(2*max_abstand_B));

	 }

  /*damit beide Molekuele dieselbe Anzahl von Atomen haben werden bei Molekuel die fehlenden Atome mit den koordinaten
      (0,0,0) initialisiert*/
    molekuel_b = (double *)realloc(molekuel_b, 3*atomzahl_a*sizeof(double));

    for(i=atomzahl_b;i<atomzahl_a;i++)
       for(j=0;j<3;j++)
           {
	     molekuel_b[3*i+j]=0.0;
	   }

     M = atomzahl_a;
     nfft_init_3d(&plan_a,64,64,64,M);
     nfft_init_3d(&plan_b,64,64,64,M);
     nfft_init_3d(&plan_a_b,64,64,64,1331);
     h_k = (complex*)fftw_malloc((plan_b.N_total)*sizeof(complex));

     for(j=0; j<plan_a.M_total; j++)
     {
         for(i=0; i<3; i++)
         {
              plan_a.x[3*j+i]=molekuel_a[3*j+i];
         }
     }

     if(plan_a.nfft_flags & PRE_ONE_PSI)
     {
          nfft_precompute_one_psi(&plan_a);
     }

     for(k=0; k<plan_a.M_total; k++)
     {
          plan_a.f[k]=gamma_A[2*k+0] + I* (gamma_A[2*k+1]);
     }

     nfft_adjoint(&plan_a);

     /*vpr_complex(plan_a.f_hat,plan_a.N_total,"adjoint nfft, vector a_hat");*/


     for(j=0;j<(this.n);j++)
     {
         this.b[j] = cexp(-this.decay)*this.R_i*csqrt(PI)/(this.p*csqrt(-this.decay))
                    *cexp((this.R_i*this.R_i*(j-(this.n)/2.0)*(j-(this.n)/2.0)*PI*PI)/(this.decay*this.p*this.p));
     }

     /* vpr_complex(this.b,this.n,"vector b_l");*/

     b_l_3d = (complex*)fftw_malloc((plan_b.N_total)*sizeof(complex));

     for(i=0;i<this.n;i++)
     {
        b_l_3d[i]=this.b[i]*this.b[i]*this.b[i];
     }

     /*  vpr_complex(b_l_3d,this.n,"vector b_l");*/


    /*Rotation des Molekuels B*/
    for(i=0;i<=360;i=i+15)
    {
        for(j=0;j<=360;j=j+15)
	{
	    for(k=0;k<=180;k=k+15)
	    {
	        alpha=i*2*PI/360;
		beta=j*2*PI/360;
		gamma=j*2*PI/360;

	        for(l=0;l<atomzahl_a;l++)
                {
                    molekuel_b_neu[3*l+0]=(cos(beta)*cos(gamma)*molekuel_b[3*l+0])
		                         +(cos(beta)*sin(gamma)*molekuel_b[3*l+1])
					 -(sin(beta)*molekuel_b[3*l+2]);
                    molekuel_b_neu[3*l+1]=((sin(alpha)*sin(beta)*cos(gamma)-(cos(alpha)*sin(gamma)))*molekuel_b[3*l+0])
		                          +((sin(alpha)*sin(beta)*sin(gamma)+(cos(alpha)*cos(gamma)))*molekuel_b[3*l+1])
					  +((sin(alpha)*cos(beta))*molekuel_b[3*l+2]);
                    molekuel_b_neu[3*l+2]=(cos(alpha)*sin(beta)*cos(gamma)+(sin(alpha)*sin(gamma)))*molekuel_b[3*l+0]
		                          +((cos(alpha)*sin(beta)*sin(gamma)-(sin(alpha)*cos(gamma)))*molekuel_b[3*l+1])
					  +((cos(alpha)*cos(beta))*molekuel_b[3*l+2]);
	        }


                for(l=0; l<plan_b.M_total; l++)
                {
                   for(r=0; r<3; r++)
                   {
                       plan_b.x[3*l+r]= molekuel_b_neu[3*l+r];
                   }
                }

                if(plan_b.nfft_flags & PRE_ONE_PSI)
                {
                    nfft_precompute_one_psi(&plan_b);
                }

                for(l=0; l<plan_b.M_total; l++)
                {
                     plan_b.f[l]=gamma_B[2*l+0] + I* (gamma_B[2*l+1]);
                }

                nfft_adjoint(&plan_b);

                /* vpr_complex(plan_b.f_hat,plan_b.N_total,"adjoint nfft, vector b_hat");*/

                h_k = (complex*)fftw_malloc((plan_b.N_total)*sizeof(complex));

                for(l=0;l<plan_b.N_total;l++)
                {
                   h_k[l]=plan_a.f_hat[l]*plan_b.f_hat[l]*b_l_3d[l]*b_l_3d[l];
                }

                l=0;
		for(s=-0.25; s<=0.25; s=s+0.05)
                {
                    for(t=-0.25;t<=0.25;t=t+0.05)
		    {
		        for(u=-0.25;u<=0.25;u=u+0.05)
			{
                           plan_a_b.x[3*l+0]=r;
			   plan_a_b.x[3*l+1]=s;
			   plan_a_b.x[3*l+2]=t;
			   l++;
			}
		     }
                }


                if(plan_a_b.nfft_flags & PRE_ONE_PSI)
                {
                     nfft_precompute_one_psi(&plan_a_b);
                }

                for(l=0; l<plan_a_b.N_total; l++)
                {
                     plan_a_b.f_hat[l]=h_k[l];
                }

                nfft_trafo(&plan_a_b);

                vpr_complex(plan_a_b.f,plan_a_b.M_total,"nfft, vector f");
          }
       }
     }

     nfft_finalize(&plan_a);

     nfft_finalize(&plan_b);

     nfft_finalize(&plan_a_b);

     free(h_k);
     free(b_l_3d);
     free(molekuel_a);
     free(molekuel_b);
     free(molekuel_b_neu);
     free(molekuel_a_extra);
     free(this.b);

}

struct octreenode* refine_tree(struct octreenode *this,double c_a[3], double l)
{
    int i;
    int j;
    double laenge;  /*Kantenlaenge des Wuefels*/
    double center[3]; /*Mittelpunkt des Wuerfels*/
    int vertices_sum;


    /*Initialisierung der Werte des uebergebenen Knotens*/
    for(i=0;i<3;i++)
    {
       this->cell_center[i] = c_a[i]; /*Zentrum des uebergebenen Wuerfels*/
    }

    this->sidelength = l; /*Kantenlaenge des uebergebenen Wuerfels*/


    /*Berechnen der Wuerfeleckpunkte*/
    for(i=0;i<8;i++)
    {
        this -> child[i] =NULL;
	if(i==0)
	{
	    this->cell_vertices[i][0]=(c_a[0]+(0.5*l));
	    this->cell_vertices[i][1]=c_a[1]+(0.5*l);
	    this->cell_vertices[i][2]=c_a[2]-(0.5*l);
	}
	 if(i==1)
	{
	   this->cell_vertices[i][0]=c_a[0]+ (0.5*l);
	   this->cell_vertices[i][1]=c_a[1]+ (0.5*l);
	   this->cell_vertices[i][2]=c_a[2]- (0.5*l);
	}
	 if(i==2)
        {
           this->cell_vertices[i][0]=c_a[0]- (0.5*l);
	   this->cell_vertices[i][1]=c_a[1]- (0.5*l);
	   this->cell_vertices[i][2]=c_a[2]- (0.5*l);
	}
	 if(i==3)
	{
	   this->cell_vertices[i][0]=c_a[0]- (0.5*l);
           this->cell_vertices[i][1]=c_a[1]+ (0.5*l);
	   this->cell_vertices[i][2]=c_a[2]- (0.5*l);
	}
	 if(i==4)
	{
	   this->cell_vertices[i][0]=c_a[0]+ (0.5*l);
           this->cell_vertices[i][1]=c_a[1]- (0.5*l);
	   this->cell_vertices[i][2]=c_a[2]+ (0.5*l);
	}
	 if(i==5)
	{
	   this->cell_vertices[i][0]=c_a[0]+ (0.5*l);
           this->cell_vertices[i][1]=c_a[1]+ (0.5*l);
	   this->cell_vertices[i][2]=c_a[2]+ (0.5*l);
	}
	 if(i==6)
	{
	   this->cell_vertices[i][0]=c_a[0]- (0.5*l);
           this->cell_vertices[i][1]=c_a[1]- (0.5*l);
	   this->cell_vertices[i][2]=c_a[2]+ (0.5*l);
	}
	 if(i==7)
	{
	   this->cell_vertices[i][0]=c_a[0]- (0.5*l);
           this->cell_vertices[i][1]=c_a[1]+ (0.5*l);
	   this->cell_vertices[i][2]=c_a[2]+ (0.5*l);
	}

    }

    /*for(i=0;i<3;i++)
    {
    printf("%lf ",this-> sidelength);
    }
    printf("\n");
    */

    /* Neuer Wert fuer die Kantenlaenge, der Kinderknoten */
      laenge=l*0.5;


        /*Rekursiver Aufruf zum Anhaengen neuer Kinder, nach vorheriger Ueberpruefung ob die Kantenlaenge der Kinder noch
	   groesser als 0.5 Angstrom ist*/

	if(laenge>2.8)
	{
	   for(i=0;i<8;i++)
	   {
              if(i==0)
	      {
	         center[0]=c_a[0]+(l*0.25);
	         center[1]=c_a[1]-(l*0.25);
	         center[2]=c_a[2]-(l*0.25);
	      }

	      if(i==1)
	      {
	         center[0]=c_a[0]+(l*0.25);
	         center[1]=c_a[1]+(l*0.25);
	         center[2]=c_a[2]-(l*0.25);
	      }

	      if(i==2)
	      {
	         center[0]=c_a[0]-(l*0.25);
	         center[1]=c_a[1]-(l*0.25);
	         center[2]=c_a[2]-(l*0.25);
	      }

	      if(i==3)
	      {
	         center[0]=c_a[0]-(l*0.25);
	         center[1]=c_a[1]+(l*0.25);
	         center[2]=c_a[2]-(l*0.25);
	      }

	      if(i==4)
	      {
	         center[0]=c_a[0]+ (l*0.25);
	         center[1]=c_a[1]- (l*0.25);
	         center[2]=c_a[2]+ (l*0.25);
	      }

	      if(i==5)
	      {
	         center[0]=c_a[0]+ (l*0.25);
	         center[1]=c_a[1]+ (l*0.25);
	         center[2]=c_a[2]+ (l*0.25);
	      }

              if(i==6)
	      {
	         center[0]=c_a[0]- (l*0.25);
	         center[1]=c_a[1]- (l*0.25);
	         center[2]=c_a[2]+(l*0.25);
	      }

	      if(i==7)
	      {
	         center[0]=c_a[0]+(l*0.25);
	         center[1]=c_a[1]+(l*0.25);
	         center[2]=c_a[2]+(l*0.25);
	      }

	      this->child[i] = (struct octreenode *) malloc(sizeof (this->child[i]));
	      this->child[i] = refine_tree(this->child[i],center,laenge);
	   }
	}
	return this;
    }

    /*traverse_tree --> Einfuegen der Atome in Core Region und Hinzufuegen neuer Atome in Skin Region*/
    void traverse_tree(struct octreenode *this, double *atomkoordinaten, double atomradius[], int molekuelzahl)
    {
	   int i,j,k,l;
	   double zentrum_atom[3]; /*Zentrum des momentanen Atoms*/
	   double abstand_mom;
	   const double vdw_radius_water = 1.8; /*van der Waals Radius von Wasser*/
	   int vertices_sum =0;
	   double vertices_mom[3];
	   int vertices_classify[8]; /*Array in dem gespeichert wird, ob die wuerfel vertices zum SAS gehoeren (sas=1)oder nicht(sas=0)*/
	   FILE *atomkoord_A_new; /*file mit neuen Atomen der Skin-Region*/
	   atomkoord_A_new= fopen("molekuel_A_new.txt","a");


	   for(i=0; i<8; i++)
	   {
	      if(this->child[i]!= NULL) /*wenn nach this->child[i] noch ein Knoten folgt d.h. der jetzige Knoten kein Blatt
	                                  ist gehe eine Ebene tiefer */
	      {
	         traverse_tree(this->child[i], atomkoordinaten, atomradius ,molekuelzahl );
	      }
	      else
	      {
	         this -> atom=0; /*Initialisierung d.h kein Atom eingefuegt*/
		 this -> sas=0; /*Initialisierung Zelle liegt außerhalb des Molekuels*/
	         for(j=0; j<molekuelzahl; j++) /*es werden alle Atome ueberprueft, ob sie im momentanen Wuerfel liegen*/
		 {
                    for(k=0; k<3; k++)
		    {
                        zentrum_atom[k]=atomkoordinaten[3*j+k];
		    }
		    abstand_mom = abstand(zentrum_atom,this-> cell_center); /*Abstand zwischen Atommittelpunkt
		                                                              Wuerfelmittelpunkt*/
		    if(abstand_mom <= atomradius[j]) /*wenn dieser Abstand kleiner als der Atomradius ist wird Atom in den
		                                       Wuerfel eingefuegt, der Atomradius dieses Atoms dem Wuerfel zugeordnet und der
						       Wuerfel zum SAS-Volume gezaehlt --> sas =1*/
		    {
                        for(k=0; k<3; k++)
			{
			    this-> atom_i[k] = zentrum_atom[k];
			    this-> radius_i =atomradius[j];
			    this -> atom=1;
			    this ->sas = 1;
			}
		    }
		 }

		 if(this->atom ==0) /*wenn jetzt noch kein kein Atom eingefuegt wurde (n->atom =0) wird ueberprueft ob der Wuerfel                                           zwischen SES- und SAS-surface liegt d.h. ob der Abstand zwischen dem Wuerfelmittelpunkt und
                                      eines Atoms kleiner ist als der Atomradius und des Radius des Wasseratoms */
		 {
		     for(j=0; j< molekuelzahl; j++)
		     {
		        for(k=0; k<3; k++)
			{
                            zentrum_atom[k]=atomkoordinaten[3*j+k];
			}
			abstand_mom = abstand(zentrum_atom,this->cell_center);
			if(abstand_mom < (atomradius[j] + vdw_radius_water))
			{
			   this -> sas=1;
			}
		     }
		 }


		 for(j=0; j<8; j++) /*Vertices der Wuerfel werden als SAS-Volume zugehoerig bzw. nicht zugehoerig 	 					klassifiziert*/
		 {
		     vertices_classify[j]=0;
		     for(k=0; k<molekuelzahl; k++)
		     {
                         for(l=0;l<3;l++)
			 {
			    zentrum_atom[l]= atomkoordinaten[3*k+l];
			    vertices_mom[l]= this-> cell_vertices[j][l];
			 }
			 abstand_mom = abstand(zentrum_atom,vertices_mom);
			 if(abstand_mom<=(atomradius[k]+ vdw_radius_water))
			 {
			     vertices_classify[j]=1; /*Wuerfeleckpunkt gehoert zum SAS*/
			 }
		     }
		     vertices_sum = vertices_sum + vertices_classify[j]; /*Summe die angibt wieviele der
		                                                                Wuerfeleckpunkte als SAS
		                                                                klassifiziert wurden, d.h. bei allen
										vertices:
										vertices_sum=8, bei keinen vertices_sum=0*/
		 }

		 if(vertices_sum<=7) /*wenn vertices_sum kleiner 7 und groesser 0 ist gehoert dieser würfel zum SAS Surface
		                       in diese
		                       Wuerfel werden jetzt neue Atome eingefuegt */
		 {
		     if(vertices_sum > 0)
		     {
		        this -> sas=2;
			this -> atom =1;
			for(k=0;k<3; k++)
			{
			   this-> atom_i[k] = this -> cell_center[k];
			}
			this -> radius_i = vdw_radius_water;

			for(j=0;j<3;j++)
			{
                           fprintf(atomkoord_A_new,"%lf ",(this ->atom_i[j]));
                        }
                        fprintf(atomkoord_A_new, "\n");

		     }
		 }
	     }
	   }
	   fclose(atomkoord_A_new);
	}

	double abstand(double vektor1[3],double vektor2[3])
	{
            double vektorabstand;
	    vektorabstand=sqrt(pow((vektor1[0]-vektor2[0]),2)
	                 +pow((vektor1[1]-vektor2[1]),2)
	                 +pow((vektor1[2]-vektor2[2]),2));
	    return vektorabstand;
        }

	void free_tree (struct octreenode *this)
	{
	    int i;
	    for(i=0;i<8;i++)
            {
                 if(this->child[i]!= NULL)
	         {
	            free_tree((this->child[i]));
	         }
	    }
	    if(this->child[0]==NULL&&this->child[1]==NULL&&this->child[2]==NULL&&this->child[3]==NULL&&this->child[4]==NULL
	       &&this->child[5]==NULL&&this->child[6]==NULL&&this->child[7]==NULL)
	    {
	       free(this);
               this=NULL;
	    }
	}


