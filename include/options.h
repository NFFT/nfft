/******************************************************************************/
/* options.h                                                                  */
/* options for the direct and fast computation of the NDFT                    */
/*                                                                            */
/* authors: D. Potts                                                          */
/*	    S. Kunis 2002                                                     */
/******************************************************************************/

/** window
 */	
#define KAISER_BESSEL
                                        /* one of KAISER_BESSEL,SINC_POWER    */
					/* B_SPLINE,GAUSSIAN,                 */

/** timing
 */
/*#define MEASURE_TIME*/
                                        /* measure time for each step         */

/** timing 
 */
#ifdef MEASURE_TIME
double elapsed_time;
#define T1 {elapsed_time=second();}
#define T2(a) {                                                                \
 elapsed_time=second()-elapsed_time;                                           \
 printf("Step %d. elapsed time: %f secs.\n",(a),elapsed_time);                 \
 /*printf("%f\t",elapsed_time);*/                                              \
 fflush(stdout);                                                               \
}
#else
#define T1
#define T2(a)
#endif

/* options.h */
