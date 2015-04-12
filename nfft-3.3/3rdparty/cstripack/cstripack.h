/*
 * Copyright (c) 2002, 2015 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#ifndef NFFT3_H
#define NFFT3_H

/* #include <f2c.h> */

int addnod_(int *nst, int *k, double *x, double *y, double *z__,
  int *list, int *lptr, int *lend, int *lnew, int *ier);

double areas_(double *v1, double *v2, double *v3);

int bdyadd_(int *kk, int *i1, int *i2, int *list, int *lptr,
  int *lend, int *lnew);

int bnodes_(int *n, int *list, int *lptr, int *lend,
  int *nodes, int *nb, int *na, int *nt);

int circum_(double *v1, double *v2, double *v3, double *c__, int *ier);

int covsph_(int *kk, int *n0, int *list, int *lptr,
  int *lend, int *lnew);

int crlist_(int *n, int *ncol, double *x, double *y, double *z__,
  int *list, int *lend, int *lptr, int *lnew, int *ltri,
  int *listc, int *nb, double *xc, double *yc, double *zc, double *rc,
  int *ier);

int delarc_(int *n, int *io1, int *io2, int * list,
  int *lptr, int *lend, int *lnew, int *ier);

int delnb_(int *n0, int *nb, int *n, int *list,
  int *lptr, int *lend, int *lnew, int *lph);

int delnod_(int *k, int *n, double *x, double *y, double *z__, int *list,
 int *lptr, int *lend, int *lnew, int *lwk, int *iwk,
 int *ier);

int edge_(int *in1, int *in2, double *x, double *y, double *z__, int *lwk,
  int *iwk, int *list, int *lptr, int *lend, int *ier);

int getnp_(double *x, double *y, double *z__, int *list, int *lptr,
  int *lend, int *l, int *npts, double *df, int *ier);

int insert_(int *k, int *lp, int *list, int *lptr,
  int *lnew);

int inside_(double *p, int *lv, double *xv, double *yv, double *zv, int *
  nv, int *listv, int *ier);

int intadd_(int *kk, int *i1, int *i2, int *i3, int *list,
  int *lptr, int *lend, int *lnew);

int intrsc_(double *p1, double *p2, double *cn, double *p, int *ier);

int jrand_(int *n, int *ix, int *iy, int *iz);

int left_(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
  double *x0, double *y0, double *z0);

int lstptr_(int *lpl, int *nb, int *list, int *lptr);

int nbcnt_(int *lpl, int *lptr);

int nearnd_(double *p, int *ist, int *n, double *x, double *y,
  double *z__, int *list, int *lptr, int *lend, double *al);

int optim_(double *x, double *y, double *z__, int *na, int *list,
  int *lptr, int *lend, int *nit, int *iwk, int *ier);

int scoord_(double *px, double *py, double *pz, double *plat, double *plon, double *pnrm);

double store_(double *x);

int swap_(int *in1, int *in2, int *io1, int *	io2,
  int *list, int *lptr, int *lend, int *lp21);

int swptst_(int *n1, int *n2, int *n3, int *n4, double *x,
  double *y, double *z__);

int trans_(int *n, double *rlat, double *rlon, double *x, double *y, double *z__);

int trfind_(int *nst, double *p, int *n, double *x, double *y, double *z__,
  int *list, int *lptr, int *lend, double *b1, double *b2,
  double *b3, int *i1, int *i2, int *i3);

int trlist_(int *n, int *list, int *lptr, int *lend,
  int *nrow, int *nt, int *ltri, int *ier);

int trlprt_(int *n, double *x, double *y, double *z__, int *iflag,
  int *nrow, int *nt, int *ltri, int *lout);

int trmesh_(int *n, double *x, double *y, double *z__, int	*list,
  int *lptr, int *lend, int *lnew, int *near__, int *next,
  double *dist, int *ier);

int trplot_(int *lun, double *pltsiz, double *elat, double *elon, double *a,
  int *n, double *x, double *y, double *z__, int *list, int *lptr,
  int *lend, char *title, int *numbr, int *ier, short title_len);

int trprnt_(int *n, double *x, double *y, double *z__, int *iflag,
  int *list, int *lptr, int *lend, int *lout);

int vrplot_(int *lun, double *pltsiz, double *elat, double *elon, double *a,
  int *n, double *x, double *y, double *z__, int *nt, int *listc,
  int *lptr, int *lend, double *xc, double *yc, double *zc, char *title,
  int *numbr, int *ier, short title_len);

#endif
