#ifndef NFFT3_H
#define NFFT3_H

/* #include <f2c.h> */

long int addnod_(long int *nst, long int *k, double *x, double *y, double *z__, 
  long int *list, long int *lptr, long int *lend, long int *lnew, long int *ier);

double areas_(double *v1, double *v2, double *v3);

long int bdyadd_(long int *kk, long int *i1, long int *i2, long int *list, long int *lptr, 
  long int *lend, long int *lnew);

long int bnodes_(long int *n, long int *list, long int *lptr, long int *lend, 
  long int *nodes, long int *nb, long int *na, long int *nt);

long int circum_(double *v1, double *v2, double *v3, double *c__, long int *ier);

long int covsph_(long int *kk, long int *n0, long int *list, long int *lptr, 
  long int *lend, long int *lnew);

long int crlist_(long int *n, long int *ncol, double *x, double *y, double *z__, 
  long int *list, long int *lend, long int *lptr, long int *lnew, long int *ltri, 
  long int *listc, long int *nb, double *xc, double *yc, double *zc, double *rc, 
  long int *ier);

long int delarc_(long int *n, long int *io1, long int *io2, long int * list, 
  long int *lptr, long int *lend, long int *lnew, long int *ier);

long int delnb_(long int *n0, long int *nb, long int *n, long int *list, 
  long int *lptr, long int *lend, long int *lnew, long int *lph);

long int delnod_(long int *k, long int *n, double *x, double *y, double *z__, long int *list, 
 long int *lptr, long int *lend, long int *lnew, long int *lwk, long int *iwk, 
 long int *ier);

long int edge_(long int *in1, long int *in2, double *x, double *y, double *z__, long int *lwk, 
  long int *iwk, long int *list, long int *lptr, long int *lend, long int *ier);

long int getnp_(double *x, double *y, double *z__, long int *list, long int *lptr, 
  long int *lend, long int *l, long int *npts, double *df, long int *ier);

long int insert_(long int *k, long int *lp, long int *list, long int *lptr, 
  long int *lnew);

long int inside_(double *p, long int *lv, double *xv, double *yv, double *zv, long int *
  nv, long int *listv, long int *ier);

long int long intadd_(long int *kk, long int *i1, long int *i2, long int *i3, long int *list, 
  long int *lptr, long int *lend, long int *lnew);

long int long intrsc_(double *p1, double *p2, double *cn, double *p, long int *ier);

long int jrand_(long int *n, long int *ix, long int *iy, long int *iz);

long int left_(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, 
  double *x0, double *y0, double *z0);

long int lstptr_(long int *lpl, long int *nb, long int *list, long int *lptr);

long int nbcnt_(long int *lpl, long int *lptr);

long int nearnd_(double *p, long int *ist, long int *n, double *x, double *y, 
  double *z__, long int *list, long int *lptr, long int *lend, double *al);

long int optim_(double *x, double *y, double *z__, long int *na, long int *list, 
  long int *lptr, long int *lend, long int *nit, long int *iwk, long int *ier);

long int scoord_(double *px, double *py, double *pz, double *plat, double *plon, double *pnrm);

double store_(double *x);

long int swap_(long int *in1, long int *in2, long int *io1, long int *	io2, 
  long int *list, long int *lptr, long int *lend, long int *lp21);

long int swptst_(long int *n1, long int *n2, long int *n3, long int *n4, double *x, 
  double *y, double *z__);

long int trans_(long int *n, double *rlat, double *rlon, double *x, double *y, double *z__);

long int trfind_(long int *nst, double *p, long int *n, double *x, double *y, double *z__, 
  long int *list, long int *lptr, long int *lend, double *b1, double *b2, 
  double *b3, long int *i1, long int *i2, long int *i3);

long int trlist_(long int *n, long int *list, long int *lptr, long int *lend, 
  long int *nrow, long int *nt, long int *ltri, long int *ier);

long int trlprt_(long int *n, double *x, double *y, double *z__, long int *iflag, 
  long int *nrow, long int *nt, long int *ltri, long int *lout);

long int trmesh_(long int *n, double *x, double *y, double *z__, long int	*list, 
  long int *lptr, long int *lend, long int *lnew, long int *near__, long int *next, 
  double *dist, long int *ier);

long int trplot_(long int *lun, double *pltsiz, double *elat, double *elon, double *a, 
  long int *n, double *x, double *y, double *z__, long int *list, long int *lptr, 
  long int *lend, char *title, long int *numbr, long int *ier, short title_len);

long int trprnt_(long int *n, double *x, double *y, double *z__, long int *iflag, 
  long int *list, long int *lptr, long int *lend, long int *lout);

long int vrplot_(long int *lun, double *pltsiz, double *elat, double *elon, double *a, 
  long int *n, double *x, double *y, double *z__, long int *nt, long int *listc, 
  long int *lptr, long int *lend, double *xc, double *yc, double *zc, char *title, 
  long int *numbr, long int *ier, short title_len);

#endif
