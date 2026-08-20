/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* Internal header for planner infrastructure. Include after "infft.h".  Never installed. */

#ifndef IPLANNER_H
#define IPLANNER_H

#include <stddef.h>
#include <stdarg.h>
#include <stdio.h>
#include <limits.h>

/* md5.c */

#if UINT_MAX >= 0xffffffffUL
typedef unsigned int md5uint;
#else
typedef unsigned long md5uint; /* >= 32 value bits per the C standard */
#endif

typedef md5uint md5sig[4];

typedef struct
{
  md5sig s;              /* running state; after md5_end: the digest words */
  unsigned char buf[64]; /* private: bytes not yet processed */
  unsigned len;          /* private: total bytes fed so far */
} md5;

void Y(md5_begin)(md5 *ctx);
void Y(md5_put_byte)(md5 *ctx, unsigned char b);
void Y(md5_put_bytes)(md5 *ctx, const void *data, size_t len);
void Y(md5_put_str)(md5 *ctx, const char *s); /* includes the '\0' */
void Y(md5_put_int)(md5 *ctx, int i);
void Y(md5_put_INT)(md5 *ctx, INT i);
void Y(md5_put_unsigned)(md5 *ctx, unsigned u);
void Y(md5_end)(md5 *ctx);

/* print.c / printers.c */
typedef struct printer_s printer;

struct printer_s {
  void (*print)(printer *p, const char *fmt, ...);
  void (*vprint)(printer *p, const char *fmt, va_list ap);
  void (*putchr)(printer *p, char c);
  void (*cleanup)(printer *p);
  int indent;
  int indent_step;
};

printer *Y(printer_create)(size_t size,
                           void (*putchr)(printer *p, char c), void (*cleanup)(printer *p));
void Y(printer_destroy)(printer *p);
printer *Y(printer_create_file)(FILE *f);
printer *Y(printer_create_str)(char *buf);
printer *Y(printer_create_cnt)(size_t *cnt);

/* scan.c / scanners.c */
typedef struct scanner_s scanner;

struct scanner_s {
  int (*scan)(scanner *sc, const char *fmt, ...);
  int (*vscan)(scanner *sc, const char *fmt, va_list ap);
  int (*getchr)(scanner *sc);
  int pushback; /* one-char pushback buffer; EOF = empty */
};

scanner *Y(scanner_create)(size_t size, int (*getchr)(scanner *sc));
void Y(scanner_destroy)(scanner *sc);
scanner *Y(scanner_create_file)(FILE *f);
scanner *Y(scanner_create_str)(const char *s);

/* primes.c */
int Y(is_prime)(INT n);
INT Y(next_prime)(INT n);

/* problem kinds */
enum {
  NFFT_PROBLEM_UNSOLVABLE,
  NFFT_PROBLEM_NFFT, // full NDFT problem
  NFFT_PROBLEM_DECONV, // problem type for forward NFFT step #1: deconvolution onto oversampled frequency grid
  NFFT_PROBLEM_CONV, // problem type for forward NFFT step #3: convolution onto nodes
  NFFT_PROBLEM_LAST
};

/* problem.c / plan.c */
typedef struct problem_s problem;
typedef struct plan_s plan;

/* problem: hashable description of work ("what") */
typedef struct
{
  int kind; // problem kind
  void (*hash)(const problem *p, md5 *ctx);
  void (*print)(const problem *p, printer *pr);
  void (*destroy)(problem *p);
} problem_adt;

struct problem_s {
  const problem_adt *adt;
};

problem *Y(problem_create)(size_t size, const problem_adt *adt);
void Y(problem_destroy)(problem *p);

/* Readiness of a plan */
enum {
  PLNR_SLEEPY = 0, // not ready
  PLNR_AWAKE_ZERO = 1, // may execute, results may be incorrect
  PLNR_AWAKE = 2 // may execute, results correct
};

typedef struct
{
  void (*apply)(const plan *ego, const problem *p); // forward transform
  void (*awake)(plan *ego, int wakefulness); // may be null
  void (*print)(const plan *ego, printer *pr);
  void (*destroy)(plan *ego); // may be null
  void (*apply_adjoint)(const plan *ego, const problem *p);  // may be null (unsupported direction)
} plan_adt;

struct plan_s {
  const plan_adt *adt;
  double pcost; // analytic estimate, may get overwritten with measured seconds
  int awake_state; // PLNR_SLEEPY | PLNR_AWAKE_ZERO | PLNR_AWAKE
  int uses_x; // whether plan reads x from problem during apply/apply_adjoint
};

plan *Y(plan_create)(size_t size, const plan_adt *adt);
void Y(plan_awake)(plan *ego, int wakefulness);
void Y(plan_destroy)(plan *ego); /* awakes to SLEEPY first */

/* timer.c */

/* Timing-loop constants */
#define PLNR_TIME_MIN_TICKS 100.0         /* accept a batch at/above this (tick mode) */
#define PLNR_TIME_MIN_SLOW_SECONDS 1.0e-3 /* accept floor in the slow-timer fallback */
#define PLNR_TIME_REPEAT 8                /* best-of batches per doubling level */
#define PLNR_TIME_LIMIT_SECONDS 2.0       /* per-candidate wall budget (budget clock) */

/* Estimate gate of the measured race: a candidate whose analytic pcost exceeds
 * PLNR_PRUNE_RATIO times the cheapest candidate's pcost is not timed. The NFFT
 * candidate set mixes complexity classes (O(N M) direct vs fast), so an ungated
 * race can spend a minute timing one apply of a direct NDFT. The ratio still
 * admits direct candidates for small problems, where the analytic models are
 * within a few x of each other. */
#define PLNR_PRUNE_RATIO 8.0

/* Measurement of pln->adt->apply(pln, p). Returns the cost in arbitrary units,
 * suitable only for comparison against other candidates measured in the same
 * process, or -1.0 when no usable clock exists. */
double Y(plan_measure_cost)(plan *pln, const problem *p);

/* solver.c */
typedef struct solver_s solver;
typedef struct planner_s planner;

typedef struct
{
  int problem_kind;
  void (*destroy)(solver *ego); /* may be NULL */
  /* return a plan with pcost set, or NULL when not applicable;
   * must be cheap and must not build node-dependent state */
  plan *(*mkplan)(const solver *ego, const problem *p, planner *pl);
} solver_adt;

struct solver_s {
  const solver_adt *adt;
  int refcnt;
};

solver *Y(solver_create)(size_t size, const solver_adt *adt); /* refcnt = 0 */
void Y(solver_use)(solver *ego);                              /* ++refcnt */
void Y(solver_destroy)(solver *ego);                          /* --refcnt; at 0: adt->destroy, free */

#define REGISTER_SOLVER(pl, s) Y(planner_register_solver)(pl, s)

/* solvtab.c */
struct solvtab_s {
  void (*reg)(planner *);
  const char *reg_nam;
};
typedef struct solvtab_s solvtab[];
void Y(solvtab_exec)(const solvtab tbl, planner *pl);
#define SOLVTAB(s)  \
  {                 \
    s, STRINGIZE(s) \
  }
#define SOLVTAB_END \
  {                 \
    0, 0            \
  }

/* planner.c */
typedef struct
{
  solver *slv;
  const char *reg_nam; /* registrar name */
  unsigned nam_hash;   /* reg_nam_hash(reg_nam), avoids strcmp on import */
  int reg_id;          /* ordinal within one registrar, 0-based */
  int next_same_kind;  /* per-kind chain, -1 terminates */
} slvdesc;

/* Impatience bounds (l, u) over the flags below; lattice order LEQ.
 * A solution with bounds (l, u) asserts: every component of the winning
 * plan is at least as impatient as l, and the search tried every solver
 * at least as impatient as u.  Solutions maintain LEQ(l, u). */
typedef struct
{
  unsigned l : 16;             /* lower impatience bound */
  unsigned info : 3;           /* PLNR_BLESSING | PLNR_H_VALID | PLNR_H_LIVE */
  unsigned timelimit_imp : 12; /* time-limit impatience */
  unsigned u : 16;             /* upper impatience bound */
  unsigned slvndx : 16;        /* winning solver index, or INFEASIBLE_SLVNDX */
} flags_t;                     /* packs into 64 bits */

/* NFFT impatience bits. Bounds convention: PLNR_ESTIMATE rides in u ONLY —
 * estimate query/memo {l = F, u = PLNR_ESTIMATE | F}, measured
 * query/bless {l = u = F}, F = forwarded user gates. ESTIMATE in l would
 * hide measured (blessed) wisdom from estimate queries. */
enum {
  PLNR_BELIEVE_PCOST = 0x0001,
  PLNR_ESTIMATE = 0x0002,
  PLNR_NO_SLOW = 0x0004,
  PLNR_NO_UGLY = 0x0008,
  PLNR_NO_DIRECT = 0x0010,     /* forbid O(N.M) direct solvers */
  PLNR_NO_NDFT_PLAIN = 0x0040, /* forbid the plain per-term 1D direct NDFT solver */
  PLNR_ALLOW_PRUNING = 0x0080,
  PLNR_NO_NDFT_BLOCKED = 0x0100, /* forbid the blocked-recurrence 1D direct NDFT solver */
  PLNR_NO_FAST_NATIVE = 0x0200   /* forbid the planner-native fast NFFT solver */
};

/* hashtable slot information */
enum {
  PLNR_BLESSING = 0x1u, /* persist this entry */
  PLNR_H_VALID = 0x2u,  /* slot has been used */
  PLNR_H_LIVE = 0x4u    /* slot holds a live entry (implies PLNR_H_VALID) */
};

/* x <= y in the impatience lattice */
#define LEQ(x, y) (((x) & (y)) == (x))
#define INFEASIBLE_SLVNDX 0xFFFFu

typedef enum {
  PLNR_FORGET_UNBLESSED,
  PLNR_FORGET_ALL
} amnesia;

typedef struct
{
  md5sig s;      /* the problem key */
  flags_t flags; /* bounds + info + winning solver */
} solution;

typedef struct
{
  solution *entries;
  unsigned size;  /* table capacity; always a prime >= 2 */
  unsigned nelem; /* live entries */
  int nrehash;    /* growth counter (observability/tests) */
} hashtab;

struct planner_s {
  slvdesc *slvdescs; /* growable array */
  unsigned nslvdesc, slvdescsiz;
  const char *cur_reg_nam; /* registration context (solvtab_exec) */
  int cur_reg_id;
  int kind_head[NFFT_PROBLEM_LAST]; /* head index per kind, -1 = empty */
  hashtab htab_blessed, htab_unblessed;
  int nthr;
  flags_t flags;
  double timelimit_seconds; /* wall-clock planning budget; -1.0 = unlimited */
};

planner *Y(planner_create)(void);
void Y(planner_destroy)(planner *pl);
void Y(planner_register_solver)(planner *pl, solver *s);
void Y(planner_forget)(planner *pl, amnesia a);
void Y(planner_export)(planner *pl, printer *p);
int Y(planner_import)(planner *pl, scanner *sc);
solution *Y(planner_hlookup)(planner *pl, const md5sig s, const flags_t *q);
void Y(planner_hinsert)(planner *pl, const md5sig s, const flags_t *f,
                        unsigned slvndx);

/* the planner's current impatience bounds */
#define PLNR_L(pl) ((pl)->flags.l)
#define PLNR_U(pl) ((pl)->flags.u)

/* Wall-clock planning timelimit. -1.0 is the "unlimited" sentinel (default).
 * The setter clamps any negative argument to -1.0 so an errant caller cannot pin the
 * planner to a zero/budgetless state. */
double Y(planner_timelimit)(const planner *pl);
void Y(planner_set_timelimit)(planner *pl, double seconds);

/* Coarse wall clock. */
double Y(planner_clock_now)(void);
double Y(planner_elapsed_seconds)(double since);

/* Returns the wisdom key. */
void Y(problem_md5)(planner *pl, const problem *p, md5sig out);

/* the estimate search: wisdom lookup, else cheapest applicable solver of
 * the problem's kind; memoises the outcome unblessed. */
plan *Y(planner_mkplan)(planner *pl, const problem *p);

/* Run every applicable solver of p's kind under the planner's current
 * bounds: store up to cap (plan, descriptor-index) pairs; return the count.
 * No wisdom lookup, no memoisation, no store mutation. */
int Y(planner_candidates)(planner *pl, const problem *p, plan **plans,
                          unsigned *slvndx, int cap);

/* Bless a winning solver for p under the planner's current bounds. */
void Y(planner_bless)(planner *pl, const problem *p, unsigned slvndx);

/* process-global planner */
planner *Y(the_planner)(void);
void Y(the_planner_destroy)(void); /* safe when absent; next call recreates */
/* creation generation of the global planner. */
unsigned Y(the_planner_generation)(void);

#define FORALL_SOLVERS(pl, s, d, what)              \
  {                                                 \
    unsigned _cnt;                                  \
    for (_cnt = 0; _cnt < (pl)->nslvdesc; ++_cnt) { \
      slvdesc *d = (pl)->slvdescs + _cnt;           \
      solver *s = d->slv;                           \
      what;                                         \
    }                                               \
  }

#define FORALL_SOLVERS_OF_KIND(kind, pl, s, d, what) \
  {                                                  \
    int _cnt = (pl)->kind_head[kind];                \
    while (_cnt >= 0) {                              \
      slvdesc *d = (pl)->slvdescs + _cnt;            \
      solver *s = d->slv;                            \
      what;                                          \
      _cnt = d->next_same_kind;                      \
    }                                                \
  }

/* tensor.c */

/* One Kronecker factor of a structured linear operator: an n_out x n_in
 * matrix acting on a strided input vector (stride is) and producing a
 * strided output vector (stride os). Note for FFTW-trained readers:
 * unlike FFTW's iodim, the input and output lengths differ in general;
 * the square case n_in == n_out recovers the iodim concept as a state,
 * not a separate type (mvdim_square, tensor_squarep). */
typedef struct
{
  INT n_in;  /* input length  (matrix columns), >= 1 */
  INT is;    /* input stride */
  INT n_out; /* output length (matrix rows), >= 1 */
  INT os;    /* output stride */
} mvdim;

/* rnk factors denote A_0 (x) ... (x) A_{rnk-1} (Kronecker product);
 * rnk == 0 is the scalar identity. Factor order is not semantically
 * significant (strides carry the layout); tensor_compress canonicalises. */
typedef struct
{
  int rnk;     /* number of factors, >= 0 */
  mvdim *dims; /* rnk entries; not used when rnk == 0 */
} tensor;

mvdim Y(mvdim_square)(INT n, INT is, INT os); /* the iodim case */

tensor *Y(tensor_create)(int rnk); /* dims allocated, uninitialised */
void Y(tensor_destroy)(tensor *t);
tensor *Y(tensor_copy)(const tensor *t);
int Y(tensor_equal)(const tensor *a, const tensor *b);
int Y(tensor_squarep)(const tensor *t);
int Y(tensor_kosherp)(const tensor *t);
INT Y(tensor_sz_in)(const tensor *t);
INT Y(tensor_sz_out)(const tensor *t);
tensor *Y(tensor_adjoint)(const tensor *t);
tensor *Y(tensor_compress)(const tensor *t);
tensor *Y(tensor_compress_contiguous)(const tensor *t);
void Y(tensor_md5)(md5 *ctx, const tensor *t);
void Y(tensor_print)(const tensor *t, printer *p);

/* NFFT kind */

typedef struct
{
  problem super;       /* kind NFFT_PROBLEM_NFFT */
  tensor *sz;          /* size of frequency tensor in caller axis order, oriented in
                          the direction of dataflow: forward N_t -> n_t; adjoint
                          is tensor_adjoint of that. */
  tensor *vecsz;       /* batching loops */
  INT M;               /* node count in time/space domain */
  int m;               /* window cutoff */
  int window;          /* NFFT_WINDOW_* ordinal */
  int sign;            /* +1 = forward, -1 = adjoint */
  unsigned fftw_flags; /* child FFTW planner flags; 0 = derive */
  R *x;                /* nodes in time/space domain. */
  C *f_hat;            /* Fourier coefficients (caller-owned) */
  C *f;                /* function values at the nodes (caller-owned) */
  int *variant;        /* per-axis NFFT_NDFT_TYPE_{I,II}; length sz->rnk, caller
                          axis order; odd axes only have one variant by definition, 
                          so normalized to TYPE_I. */
  int x_owned;         /* 1: x is this problem's private copy;
                          0: x is borrowed from a parent problem. */
} problem_nfft;

problem *Y(mkproblem_nfft)(int d, const INT *N, const int *variant,
                           const INT *n, INT M, int m, int window, int sign,
                           unsigned fftw_flags, R *x, int copy_x, C *f_hat, C *f);

/* direction-aware accessors */
INT Y(problem_nfft_N)(const problem *p, int t);
INT Y(problem_nfft_n)(const problem *p, int t);
INT Y(problem_nfft_Ntot)(const problem *p);
INT Y(problem_nfft_ntot)(const problem *p);

/* per-axis NDFT variant (NFFT_NDFT_TYPE_{I,II}); t in [0, sz->rnk). */
int Y(problem_nfft_variant)(const problem *p, int t);

/* 1 if any axis has N_t == 1. Unit axes are elided at construction, so a
 * solver that sees one was reached with compression bypassed; every rank >= 1
 * solver declines. A rank-0 problem has no axes and is served by the ungated
 * base case. */
int Y(problem_nfft_has_unit_axis)(const problem *p);

/* x-restore byte-identity guard (kernel/nfft/xcheck.c): an md5 checksum of an
 * rnk*M-length real array */
void Y(nfft_x_md5)(const R *x, INT dM, md5sig out);
int Y(nfft_x_verify)(const R *x, INT dM, const md5sig ref);

/* NFFT solver roster. */
void Y(nfft_solvers_register)(planner *pl);
void Y(nfft_ensure_registered)(void);

/* Zero-dimensional direct NDFT aka the constant function.
* forward = broadcast, adjoint = reduce. */
void Y(nfft_solver_rnk0_register)(planner *pl);

/* One-dimensional direct NDFT */
void Y(nfft_solver_ndft_1d_register)(planner *pl);

/* One-dimensional blocked direct NDFT */
void Y(nfft_solver_ndft_1d_blocked_register)(planner *pl);

/* Multi-dimensional direct NDFT */
void Y(nfft_solver_ndft_nd_register)(planner *pl);

/* NFFT, decomposes problem into DECONV, FFTW, and CONV children. */
void Y(nfft_solver_fast_native_register)(planner *pl);

/* Solver plans. */
struct Y(plan_ng_s);
void Y(plan_ng_print)(struct Y(plan_ng_s) * p, printer *pr);

/* The compile-time-selected window. */
int Y(get_window_id)(void);

/* Runtime window evaluation: k = actual frequency index; x = signed grid distance. */
R Y(window_phi_hut)(int window, INT n, INT N, int m, INT k);
R Y(window_phi)(int window, INT n, INT N, int m, R x);

/*   - phi_hut_apply: DECONV Step A, regular frequency band k0..k0+count-1.
 *   - phi_precompute: CONV Step C, per-axis PRE_PSI psi -- owns c=floor(n*x)
 *     centering + the 2m+2 tap loop, writing the strided psi layout directly. */
void Y(window_phi_hut_apply)(int window, INT n, INT N, int m, INT k0,
                             R *out, INT count);
void Y(window_phi_precompute)(int window, INT n, INT N, int m,
                              const R *x, INT x_stride, INT num_nodes,
                              R *out, INT out_stride);

/* NFFT deconvolution.  Maps the input frequency tensor onto the oversampled grid, 
 * dividing by the window's Fourier coefficients. Step is node-independent. */
typedef struct
{
  problem super; /* kind NFFT_PROBLEM_DECONV */
  tensor *sz;    /* size of frequency tensor in caller axis order, oriented in
                          the direction of dataflow: forward N_t -> n_t; adjoint
                          is tensor_adjoint of that. */
  tensor *vecsz; /* batching loops */
  int m;         /* window cutoff */
  int window;    /* NFFT_WINDOW_* ordinal */
  int sign;      /* +1 = deconvolve+zero-pad f_hat->g, -1 = adjoint g->f_hat */
  int *variant;  /* per-axis NFFT_NDFT_TYPE_{I,II}; length sz->rnk, caller
                    axis order; odd axes only have one variant by definition, 
                    so normalized to TYPE_I. */
  C *f_hat;      /* borrowed in (parent's f_hat) */
  C *g;          /* borrowed out */
} problem_deconv;

problem *Y(mkproblem_deconv)(int d, const INT *N, const int *variant,
                             const INT *n, int m, int window, int sign,
                             C *f_hat, C *g);
INT Y(problem_deconv_N)(const problem *p, int t);
INT Y(problem_deconv_n)(const problem *p, int t);
INT Y(problem_deconv_Ntot)(const problem *p);
INT Y(problem_deconv_ntot)(const problem *p);
int Y(problem_deconv_variant)(const problem *p, int t);

/* DECONV solver roster */
void Y(deconv_solvers_register)(planner *pl);
void Y(deconv_ensure_registered)(void);

/* NFFT convolution Maps the oversampled grid to the nonequispaced node samples 
 * via the window convolution. */
typedef struct
{
  problem super; /* kind NFFT_PROBLEM_CONV */
  tensor *sz;    /* square grid n_t -> n_t */
  tensor *vecsz; /* batching loops */
  INT *N;        /* per-axis bandwidth; geometry, needed for the KB shape
                  * parameter b_t = pi(2 - N_t/n_t) */
  INT M;         /* node count */
  int m;         /* window cutoff */
  int window;    /* NFFT_WINDOW_* ordinal */
  int sign;      /* +1 = g->f, -1 = adjoint f->g */
  R *x;          /* caller-owned nodes */
  C *g;          /* borrowed in */
  C *f;          /* borrowed out */
} problem_conv;

problem *Y(mkproblem_conv)(int d, const INT *n, const INT *N, INT M, int m,
                           int window, int sign, R *x, C *g, C *f);
INT Y(problem_conv_n)(const problem *p, int t);
INT Y(problem_conv_N)(const problem *p, int t);
INT Y(problem_conv_ntot)(const problem *p);

/* CONV solver roster */
void Y(conv_solvers_register)(planner *pl);
void Y(conv_ensure_registered)(void);

#endif /* IPLANNER_H */
