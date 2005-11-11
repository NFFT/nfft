/*****************************************************************************/
/* options.h                                                                 */
/* options for the direct and fast computation of the NDFT                   */
/*                                                                           */
/* authors: D. Potts                                                         */
/*	    S. Kunis 2002                                                    */
/*****************************************************************************/

/** window
 */	
#define GAUSSIAN
                                        /* one of KAISER_BESSEL,SINC_POWER   */
					/* B_SPLINE,GAUSSIAN,                */

/** timing
 */
//#define MEASURE_TIME
                                        /* measure time for each step        */

/** timing 
 */
#ifdef MEASURE_TIME
 int MEASURE_TIME_r;
 double MEASURE_TIME_tt;

#define TIC(a)                                                                \
  ths->MEASURE_TIME_t[(a)]=0;                                                 \
  MEASURE_TIME_r=0;                                                           \
  while(ths->MEASURE_TIME_t[(a)]<0.01)                                        \
    {                                                                         \
      MEASURE_TIME_r++;                                                       \
      MEASURE_TIME_tt=second();                                               \

/* THE MEASURED FUNCTION IS CALLED REPEATEDLY */

#define TOC(a)                                                                \
      MEASURE_TIME_tt=second()-MEASURE_TIME_tt;                               \
      ths->MEASURE_TIME_t[(a)]+=MEASURE_TIME_tt;                              \
    }                                                                         \
  ths->MEASURE_TIME_t[(a)]/=MEASURE_TIME_r;                                   \

#else
#define TIC(a)
#define TOC(a)
#endif

/* options.h */
