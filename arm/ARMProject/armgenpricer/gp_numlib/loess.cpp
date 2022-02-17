/*!
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *	\file loess.cpp
 *
 *  \brief 
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#include "gpnumlib\loess.h"
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <string>
#include "expt.h"

CC_BEGIN_NAMESPACE( ARM )

#define Calloc(n,t)	(t *)calloc((unsigned)(n),sizeof(t))
#define Free(p)		free((char *)(p))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define abs(x) ((x) >= 0 ? (x) : -(x))

#define	GAUSSIAN	1
#define SYMMETRIC	0

// Intermediate Functions

void loess_raw(
double *y, 
double *x, 
double *weights, 
double *robust, 
long *d, 
long *n, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
double *cell, 
char **surf_stat, 
double *surface, 
long *parameter, 
long *a, 
double *xi, 
double *vert, 
double *vval, 
double *diagonal, 
double *trL, 
double *one_delta, 
double *two_delta, 
long *setLf);

void loess_workspace(
long *d, 
long *n, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
long *setLf);

void loess_prune(
long *parameter, 
long *a, 
double *xi, 
double *vert, 
double *vval);

void loess_grow(
long *parameter, 
long *a, 
double *xi, 
double *vert, 
double *vval);

void loess_dfitse(
double *y,
double *x, 
double *x_evaluate, 
double *weights, 
double *robust, 
long *family, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
long *d, 
long *n, 
long *m, 
double *fit, 
double *L);

void loess_dfit(
double *y, 
double *x, 
double *x_evaluate, 
double *weights, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
long *d, 
long *n, 
long *m, 
double *fit);

void loess_ifit(
double *parameter, 
double *a,
double *xi, 
double *vert, 
double *vval, 
double *m, 
double *x_evaluate, 
double *fit);

void loess_ise(
double *y, 
double *x, 
double *x_evaluate, 
double *weights, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
double *cell, 
long *d, 
long *n, 
long *m, 
double *fit, 
double *L);


int lowesa(double *trl, long *n, long *d, long 
	*tau, long *nsing, double *delta1, double *delta2);

int lowesb(double *xx, double *yy, double *ww, 
	double *diagl, long *infl, long *iv, long *liv, long *
	lv, double *wv);

int lowesc(long *n, double *l, double *ll, 
	double *trl, double *delta1, double *delta2);

int lowesd(long *versio, long *iv, long *liv, 
	long *lv, double *v, long *d, long *n, double *f, 
	long *ideg, long *nvmax, long *setlf);

int lowese(long *iv, long *liv, long *lv, 
	double *wv, long *m, double *z, double *s);

int lowesf(double *xx, double *yy, double *ww, 
	long *iv, long *liv, long *lv, double *wv, long *m, 
	double *z, double *l, long *ihat, double *s);

int lowesl(long *iv, long *liv, long *lv, 
	double *wv, long *m, double *z, double *l);

int lowesp(long *n, double *y, double *yhat, double *pwgts, double *rwgts, long *pi, double *ytilde);


int lowesw(double *res, long *n, double *rw, long *pi);


int ehg169_(long *d, long *vc, long *nc, long *
	ncmax, long *nv, long *nvmax, double *v, long *a, 
	double *xi, long *c, long *hi, long *lo);

int ehg196(long *tau, long *d, double *f, 
	double *trl);

int ehg139_(double *v, long *nvmax, long *nv, 
	long *n, long *d, long *nf, double *f, double *x, 
	long *pi, long *psi, double *y, double *rw, double *
	trl, long *kernel, long *k, double *dist, long *phi, 
	double *eta, double *b, long *od, double *w, 
	double *diagl, double *vval2, long *ncmax, long *vc, 
	long *a, double *xi, long *lo, long *hi, long *c, 
	long *vhit, double *rcond, long *sing, long *dd, long 
	*tdeg, long *cdeg, long *lq, double *lf, long *setlf, 
	double *s);

typedef struct
{	long cierr;
	long ciunit;
	long ciend;
	char *cifmt;
	long cirec;
} cilist;

// Fortran Function

double gamma(double xx);
double pf(double q, double df1, double df2);
double qt(double p, double df);
double ibeta(double x, double a, double b);
double invigauss_quick(double p);
double invibeta(double p, double a, double b);
double invibeta_quick(double p, double a, double b);
double fmin(double a, double b);
double fmax(double a, double b);
void Recover(char* a, int *b);
void Warning(char *a, int *b);

double ddot_(long *n, double *dx, long *incx, double *dy, long *incy);
double d1mach_(long *i);
int dqrsl_(double *x, long *ldx, long *n, long *k, double *qraux, double *y, double *qy, double *qty, double *b, double *rsd, double *xb, long *job, long *info);
int dcopy_(long *n, double *dx, long *incx, double *dy, long *incy);
int daxpy_(long *n, double *da, double *dx, long *incx, double *dy, long *incy);
int dsvdc_(double *x, int *ldx, int *n, int *p, double *s, double *e, double *u, int *ldu, double *v, int *ldv, double *work, int *job, int *info);
int drot_(long *n, double *dx, long *incx, double *dy, long *incy, double *c, double *s);
int drotg_(double *da, double *db, double *c, double *s);
int dswap_(long *n, double *dx, long *incx, double *dy, long *incy);
int dscal_(long *n, double *da, double *dx, long *incx);
double dnrm2_(long *, double* x, long* incx);
extern double c_b44;

extern "C" double d_sign(double*, double*);
extern "C" double pow_dd(double *,double *);
extern "C" long pow_ii(long *,long *);
extern "C" int s_stop(char *,long);
extern "C" long e_wsle(void);
extern "C" long s_wsle(cilist *);
extern "C" long do_lio(long *,long *,char *,long);
extern "C" long s_wsfe(cilist *);
extern "C" long e_wsfe(void);
extern "C" long e_rsle(void);
extern "C" long s_rsle(cilist *);
extern "C" long do_fio(long *,char *,long);
extern "C" double sqrt(double);
extern "C" double exp(double);
extern "C" double log(double);
extern "C" double log10(double);
extern "C" double fabs(double);
extern "C" double ceil(double);
extern "C" double floor(double);
extern "C" double pow(double,double);

extern "C" void* calloc(size_t, size_t);
extern "C" void free(void *);
extern "C" void* malloc(size_t);
extern "C" void qsort(void *, size_t, size_t, int (__cdecl *)(const void *,const void *));
extern "C" void exit(int);

static long	*iv, liv, lv, tau;
static double	*v;

// Loess Functions

void loess_raw(
double *y, 
double *x, 
double *weights, 
double *robust, 
long *d, 
long *n, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
double *cell, 
char **surf_stat, 
double *surface, 
long *parameter, 
long *a, 
double *xi, 
double *vert, 
double *vval, 
double *diagonal, 
double *trL, 
double *one_delta, 
double *two_delta, 
long *setLf)
{
	long	zerol = 0, one = 1, two = 2, nsing, i, k;
	double	*hat_matrix, *LL, zero = 0.0;


	*trL = 0;
	loess_workspace(d, n, span, degree, nonparametric, drop_square, 
		sum_drop_sqr, setLf);
        v[1] = *cell;
	if(!strcmp(*surf_stat, "interpolate/none")) {
		lowesb(x, y, robust, &zero, &zerol, iv, &liv, &lv, v);
		lowese(iv, &liv, &lv, v, n, x, surface);
		loess_prune(parameter, a, xi, vert, vval);
	}			
	else if (!strcmp(*surf_stat, "direct/none")) {
		lowesf(x, y, robust, iv, &liv, &lv, v, n, x,
			&zero, &zerol, surface);
	}
	else if (!strcmp(*surf_stat, "interpolate/1.approx")) {
		lowesb(x, y, weights, diagonal, &one, iv, &liv, &lv, v);
		lowese(iv, &liv, &lv, v, n, x, surface);
		nsing = iv[29];
		for(i = 0; i < (*n); i++) *trL = *trL + diagonal[i];
		lowesa(trL, n, d, &tau, &nsing, one_delta, two_delta);
		loess_prune(parameter, a, xi, vert, vval);
	}
        else if (!strcmp(*surf_stat, "interpolate/2.approx")) {
		lowesb(x, y, robust, &zero, &zerol, iv, &liv, &lv, v);
		lowese(iv, &liv, &lv, v, n, x, surface);
		nsing = iv[29];
		ehg196(&tau, d, span, trL);
		lowesa(trL, n, d, &tau, &nsing, one_delta, two_delta);
		loess_prune(parameter, a, xi, vert, vval);
	}
	else if (!strcmp(*surf_stat, "direct/approximate")) {
		lowesf(x, y, weights, iv, &liv, &lv, v, n, x,
			diagonal, &one, surface);
		nsing = iv[29];
		for(i = 0; i < (*n); i++) *trL = *trL + diagonal[i];
		lowesa(trL, n, d, &tau, &nsing, one_delta, two_delta);
	}
	else if (!strcmp(*surf_stat, "interpolate/exact")) {
		hat_matrix = Calloc((*n)*(*n), double);
		LL = Calloc((*n)*(*n), double);
		lowesb(x, y, weights, diagonal, &one, iv, &liv, &lv, v);
		lowesl(iv, &liv, &lv, v, n, x, hat_matrix);
		lowesc(n, hat_matrix, LL, trL, one_delta, two_delta);
		lowese(iv, &liv, &lv, v, n, x, surface);
		loess_prune(parameter, a, xi, vert, vval);
		Free(hat_matrix);
		Free(LL);
	}
	else if (!strcmp(*surf_stat, "direct/exact")) {
		hat_matrix = Calloc((*n)*(*n), double);
		LL = Calloc((*n)*(*n), double);
		lowesf(x, y, weights, iv, &liv, &lv, v, n, x,
			hat_matrix, &two, surface);
		lowesc(n, hat_matrix, LL, trL, one_delta, two_delta);
                k = (*n) + 1;
		for(i = 0; i < (*n); i++)
			diagonal[i] = hat_matrix[i * k];
		Free(hat_matrix);
		Free(LL);
	}
	loess_free();
}

void loess_dfit(
double *y, 
double *x, 
double *x_evaluate, 
double *weights, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
long *d, 
long *n, 
long *m, 
double *fit)
{
	long	zerol = 0, one = 1;
	double zero = 0.0;
	
        loess_workspace(d, n, span, degree, nonparametric, drop_square,
                sum_drop_sqr, &zerol);
	lowesf(x, y, weights, iv, &liv, &lv, v, m, x_evaluate,
			&zero, &zerol, fit);
	loess_free();
}

void loess_dfitse(
double *y,
double *x, 
double *x_evaluate, 
double *weights, 
double *robust, 
long *family, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
long *d, 
long *n, 
long *m, 
double *fit, 
double *L)
{
	long	zerol = 0, one = 1, two = 2;
	double zero = 0.0;
	
        loess_workspace(d, n, span, degree, nonparametric, drop_square,
                sum_drop_sqr, &zerol);
	if(*family == GAUSSIAN)
		lowesf(x, y, weights, iv, &liv, &lv, v, m, 
				x_evaluate, L, &two, fit);
	else if(*family == SYMMETRIC)
	{
		lowesf(x, y, weights, iv, &liv, &lv, v, m,
				x_evaluate, L, &two, fit);
		lowesf(x, y, robust, iv, &liv, &lv, v, m,
				x_evaluate, &zero, &zerol, fit);
	}	
	loess_free();
}
void loess_ifit(
double *parameter, 
double *a,
double *xi, 
double *vert, 
double *vval, 
double *m, 
double *x_evaluate, 
double *fit)
{
	long* ml = (long*) m;
	long *parameterl = (long*) parameter;
	long *la = (long*) a;
	loess_grow(parameterl, la, xi, vert, vval);
	lowese(iv, &liv, &lv, v, ml, x_evaluate, fit);
	loess_free();
}

void loess_ise(
double *y, 
double *x, 
double *x_evaluate, 
double *weights, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
double *cell, 
long *d, 
long *n, 
long *m, 
double *fit, 
double *L)
{
	long	zerol = 0, one = 1;
	double zero = 0.0;
	
    loess_workspace(d, n, span, degree, nonparametric, drop_square,
                sum_drop_sqr, &one);
	v[1] = *cell;
	lowesb(x, y, weights, &zero, &zerol, iv, &liv, &lv, v);
	lowesl(iv, &liv, &lv, v, m, x_evaluate, L);
	loess_free();
}

void loess_workspace(
long *d, 
long *n, 
double *span, 
long *degree, 
long *nonparametric, 
long *drop_square, 
long *sum_drop_sqr, 
long *setLf)
{
	long	D, N, tau0, nvmax, nf, version = 106, i;

	D = *d;
	N = *n;
	nvmax = max(200, N);
        nf = (long) min(N, floor(N * (*span)));
        tau0 = (long) (((*degree) > 1) ? ((D + 2) * (D + 1) * 0.5) : (D + 1));
        tau = tau0 - (*sum_drop_sqr);
        lv = 50 + (3 * D + 3) * nvmax + N + (tau0 + 2) * nf;
	liv = 50 + ((long)pow((double)2, (double)D) + 4) * nvmax + 2 * N;
	if(*setLf) {
		lv = lv + (D + 1) * nf * nvmax;
		liv = liv + nf * nvmax;	
	}
        iv = Calloc(liv, long);
        v = Calloc(lv, double);

        lowesd(&version, iv, &liv, &lv, v, d, n, span, degree, 
			&nvmax, setLf);
        iv[32] = *nonparametric;
        for(i = 0; i < D; i++)
                iv[i + 40] = drop_square[i];
}

void loess_prune(
long *parameter, 
long *a, 
double *xi, 
double *vert, 
double *vval)
{
	long	d, vc, a1, v1, xi1, vv1, nc, nv, nvmax, i, k;
	
	d = iv[1];
	vc = iv[3] - 1;
	nc = iv[4];
	nv = iv[5];
	a1 = iv[6] - 1;
	v1 = iv[10] - 1;
	xi1 = iv[11] - 1;
	vv1 = iv[12] - 1;
	nvmax = iv[13];

	for(i = 0; i < 5; i++)
		parameter[i] = iv[i + 1];
	parameter[5] = iv[21] - 1;
	parameter[6] = iv[14] - 1;

	for(i = 0; i < d; i++){
		k = nvmax * i;
		vert[i] = v[v1 + k];
		vert[i + d] = v[v1 + vc + k];
	}
	for(i = 0; i < nc; i++) {
		xi[i] = v[xi1 + i];
		a[i] = iv[a1 + i];
	}
	k = (d + 1) * nv;
	for(i = 0; i < k; i++)
		vval[i] = v[vv1 + i];
}

void loess_grow(
long *parameter, 
long *a, 
double *xi, 
double *vert, 
double *vval)
{
	long	d, vc, nc, nv, a1, v1, xi1, vv1, i, k;

	d = parameter[0];
	vc = parameter[2];
	nc = parameter[3];
	nv = parameter[4];
	liv = parameter[5];
	lv = parameter[6];
	iv = Calloc(liv, long);
	v = Calloc(lv, double);

	iv[1] = d;
	iv[2] = parameter[1];
	iv[3] = vc;
	iv[5] = iv[13] = nv;
	iv[4] = iv[16] = nc;
	iv[6] = 50;
	iv[7] = iv[6] + nc;
	iv[8] = iv[7] + vc * nc;
	iv[9] = iv[8] + nc;
	iv[10] = 50;
	iv[12] = iv[10] + nv * d;
	iv[11] = iv[12] + (d + 1) * nv;
	iv[27] = 173;

	v1 = iv[10] - 1;
	xi1 = iv[11] - 1;
	a1 = iv[6] - 1;
	vv1 = iv[12] - 1;
	
        for(i = 0; i < d; i++) {
		k = nv * i;
		v[v1 + k] = vert[i];
		v[v1 + vc - 1 + k] = vert[i + d];
	}
        for(i = 0; i < nc; i++) {
                v[xi1 + i] = xi[i];
                iv[a1 + i] = a[i];
        }
	k = (d + 1) * nv;
	for(i = 0; i < k; i++)
		v[vv1 + i] = vval[i];

	ehg169_(&d, &vc, &nc, &nc, &nv, &nv, v+v1, iv+a1,
			v+xi1, iv+iv[7]-1, iv+iv[8]-1, iv+iv[9]-1);
}

void loess_free()
{
        Free(v);
        Free(iv);
}

// Predict Functions

void pred_(
double *y, 
double  *x_, 
double  *new_x, 
long *size_info, 
double  *s, 
double  *weights, 
double  *robust, 
double  *span, 
long *degree, 
long *normalize, 
long *parametric, 
long *drop_square, 
char **surface, 
double *cell, 
char **family, 
long *parameter, 
long *a, 
double  *xi, 
double  *vert, 
double  *vval, 
double  *divisor, 
long *se, 
double  *fit, 
double  *se_fit)
{     
        double  *x, *x_tmp, *x_evaluate, *L, new_cell, tmp, *fit_tmp, 
	        *temp;
        long    N, D, M, sum_drop_sqr = 0, sum_parametric = 0,
	        nonparametric = 0, *order_parametric, *order_drop_sqr;
	int     i, j, k, p, comp();

        D = size_info[0];
        N = size_info[1];
	M = size_info[2];

	x = (double *) malloc(N * D * sizeof(double));
	x_tmp = (double *) malloc(N * D * sizeof(double));
	x_evaluate = (double *) malloc(M * D * sizeof(double));
	L = (double *) malloc(N * M * sizeof(double));
        order_parametric = (long *) malloc(D * sizeof(long));
        order_drop_sqr = (long *) malloc(D * sizeof(long));
	temp = (double *) malloc(N * D * sizeof(double));

	for(i = 0; i < (N * D); i++)
		x_tmp[i] = x_[i];
	for(i = 0; i < D; i++) {
		k = i * M;
		for(j = 0; j < M; j++) {
			p = k + j;
			new_x[p] = new_x[p] / divisor[i];
		}
	}
	if(!strcmp(*surface, "direct") || se) {
		for(i = 0; i < D; i++) {
			k = i * N;
			for(j = 0; j < N; j++) {
                                p = k + j;
                                x_tmp[p] = x_[p] / divisor[i];
                        }
		}
	}
	j = D - 1;
	for(i = 0; i < D; i++) {
	        sum_drop_sqr = sum_drop_sqr + drop_square[i];
	        sum_parametric = sum_parametric + parametric[i];
	        if(parametric[i])
                        order_parametric[j--] = i;
		else
	                order_parametric[nonparametric++] = i;
	}
        for(i = 0; i < D; i++) {
                order_drop_sqr[i] = 2 - drop_square[order_parametric[i]];
		k = i * M;
		p = order_parametric[i] * M;
	        for(j = 0; j < M; j++)
			x_evaluate[k + j] = new_x[p + j];
		k = i * N;
		p = order_parametric[i] * N;
	        for(j = 0; j < N; j++)
		        x[k + j] = x_tmp[p + j];
        }
	for(i = 0; i < N; i++)
		robust[i] = weights[i] * robust[i];

	if(!strcmp(*surface, "direct")) {
	        if(*se) {
				long test = !strcmp(*family, "gaussian");
		        loess_dfitse(y, x, x_evaluate, weights, robust,
				&test, span, degree,
                                &nonparametric, order_drop_sqr, &sum_drop_sqr,
                                &D, &N, &M, fit, L);
                }
	        else {
		        loess_dfit(y, x, x_evaluate, robust, span, degree,
                                &nonparametric, order_drop_sqr, &sum_drop_sqr,
				&D, &N, &M, fit);
                }
        }
	else {
			double *dM = (double*) &M;
			double *da = (double*) a;
			double *dparameter = (double*) parameter;
	        loess_ifit(dparameter, da, xi, vert, vval, dM, x_evaluate, fit);
	        if(*se) {
                        new_cell = (*span) * (*cell);
           	        fit_tmp = (double *) malloc(M * sizeof(double));
		        loess_ise(y, x, x_evaluate, weights, span, degree,
				&nonparametric, order_drop_sqr, &sum_drop_sqr,
				&new_cell, &D, &N, &M, fit_tmp, L);
			free(fit_tmp);
                }
        }
	if(*se) {
	        for(i = 0; i < N; i++) {
		        k = i * M;
	                for(j = 0; j < M; j++) {
			        p = k + j;
			        L[p] = L[p] / weights[i];
			        L[p] = L[p] * L[p];
			}
		}
		for(i = 0; i < M; i++) {
		        tmp = 0;
			for(j = 0; j < N; j++)
			        tmp = tmp + L[i + j * M];
			se_fit[i] = (*s) * sqrt(tmp);
		}
	}
	free(x);
	free(x_tmp);
	free(x_evaluate);
	free(L);
	free(order_parametric);
	free(order_drop_sqr);
	free(temp);
}

void predict(
double *eval, 
long m, 
struct loess_struct *lo, 
struct pred_struct *pre, 
long se)
{
	long	size_info[3];

        pre->fit = (double *) malloc(m * sizeof(double));
        pre->se_fit = (double *) malloc(m * sizeof(double));
	pre->residual_scale = lo->out.s;
	pre->df = (lo->out.one_delta * lo->out.one_delta) / lo->out.two_delta;

	size_info[0] = lo->in.p;
	size_info[1] = lo->in.n;
	size_info[2] = m;
	
	pred_(lo->in.y, lo->in.x, eval, size_info, &lo->out.s, 
		lo->in.weights, 
		lo->out.robust,
		&lo->model.span,
		&lo->model.degree,
		&lo->model.normalize,
		lo->model.parametric,
		lo->model.drop_square,
		&lo->control.surface,
		&lo->control.cell,
		&lo->model.family,
		lo->kd_tree.parameter,
		lo->kd_tree.a,
		lo->kd_tree.xi,
		lo->kd_tree.vert,
		lo->kd_tree.vval,
		lo->out.divisor,
		&se,
		pre->fit,
		pre->se_fit);
}

void pred_free_mem(struct pred_struct	*pre)
{
	free(pre->fit);
	free(pre->se_fit);
}

void ehg183(char *s, int *i, int *n,int *inc)

{
  char mess[4000], num[20];
  int j;
  strcpy(mess,s);
  for (j=0; j<*n; j++) {
    sprintf(num," %d",i[j * *inc]);
    strcat(mess,num);
  }
  strcat(mess,"\n");
  Warning(mess,NULL);
}

void ehg184(char *s, double *x, int *n, int *inc)
{
  char mess[4000], num[30];
  int j;
  strcpy(mess,s);
  for (j=0; j<*n; j++) {
    sprintf(num," %.5g",x[j * *inc]);
    strcat(mess,num);
  }
  strcat(mess,"\n");
  Warning(mess,NULL);
}


static  char    *surf_stat;

void loess_setup(double *x, double *y, long n, long p, struct  loess_struct* lo)
{
	int	i, max_kd;

	max_kd = n > 200 ? n : 200;

	lo->in.y = (double *) malloc(n * sizeof(double));
        lo->in.x = (double *) malloc(n * p * sizeof(double));
	lo->in.weights = (double *) malloc(n * sizeof(double));
	for(i = 0; i < (n * p); i++)
	        lo->in.x[i] = x[i];
	for(i = 0; i < n; i++) {
	        lo->in.y[i] = y[i];
		lo->in.weights[i] = 1;
	}
	lo->in.n = n;
	lo->in.p = p;
        lo->model.span = 0.75;
	lo->model.degree = 2;
	lo->model.normalize = true;
	for(i = 0; i < 8; i++)
	        lo->model.parametric[i] = lo->model.drop_square[i] = false;
	lo->model.family = "gaussian";
        lo->control.surface = "interpolate";
        lo->control.statistics = "approximate";
	lo->control.cell = 0.2;
	lo->control.trace_hat = "wait.to.decide";
 	lo->control.iterations = 4;

	lo->out.fitted_values = (double *) malloc(n * sizeof(double));
	lo->out.fitted_residuals = (double *) malloc(n * sizeof(double));
	lo->out.pseudovalues = (double *) malloc(n * sizeof(double));
	lo->out.diagonal = (double *) malloc(n * sizeof(double));
	lo->out.robust = (double *) malloc(n * sizeof(double));
	lo->out.divisor = (double *) malloc(p * sizeof(double));

	lo->kd_tree.parameter = (long *) malloc(7 * sizeof(long));
	lo->kd_tree.a = (long *) malloc(max_kd * sizeof(long));
	lo->kd_tree.xi = (double *) malloc(max_kd * sizeof(double));
	lo->kd_tree.vert = (double *) malloc(p * 2 * sizeof(double));
	lo->kd_tree.vval = (double *) malloc((p + 1) * max_kd * sizeof(double));
}	

void condition(char **surface, char *new_stat, char **trace_hat_in)
{
	if(!strcmp(*surface, "interpolate")) {
		if(!strcmp(new_stat, "none"))
			surf_stat = "interpolate/none";
		else if(!strcmp(new_stat, "exact"))
			surf_stat = "interpolate/exact";
		else if(!strcmp(new_stat, "approximate"))
		{
			if(!strcmp(*trace_hat_in, "approximate"))
				surf_stat = "interpolate/2.approx";
			else if(!strcmp(*trace_hat_in, "exact"))
				surf_stat = "interpolate/1.approx";
		}
	}
	else if(!strcmp(*surface, "direct")) {
		if(!strcmp(new_stat, "none"))
			surf_stat = "direct/none";
		else if(!strcmp(new_stat, "exact"))
			surf_stat = "direct/exact";
		else if(!strcmp(new_stat, "approximate"))
			surf_stat = "direct/approximate";
	}
}

int comp(double *d1, double *d2)
{
        if(*d1 < *d2)
                return(-1);
        else if(*d1 == *d2)
                return(0);
        else
                return(1);
}

int vcomp(const void *v1, const void *v2)
{
	return comp((double*) v1, (double*) v2);
}

void loess_(
double *y, 
double *x_, 
long *size_info, 
double *weights, 
double *span, 
long *degree, 
long *parametric, 
long *drop_square,
long *normalize, 
char **statistics, 
char **surface, 
double *cell, 
char **trace_hat_in, 
long *iterations,
double *fitted_values, 
double *fitted_residuals, 
double *enp, 
double *s, 
double *one_delta, 
double *two_delta, 
double *pseudovalues, 
double *trace_hat_out, 
double *diagonal, 
double *robust, 
double *divisor, 
long *parameter, 
long *a, 
double *xi, 
double *vert, 
double *vval)
{
	double	*x, *x_tmp, new_cell, trL, delta1, delta2, sum_squares = 0, 
		*pseudo_resid, *temp, *xi_tmp, *vert_tmp, *vval_tmp, 
		*diag_tmp, trL_tmp = 0, d1_tmp = 0, d2_tmp = 0, sum, mean;
	long	i, j, k, p, N, D, sum_drop_sqr = 0, sum_parametric = 0, 
		setLf,	nonparametric = 0, *order_parametric,
		*order_drop_sqr, zero = 0, max_kd, *a_tmp, *param_tmp;
	int     cut;
	char	*new_stat;

	D = size_info[0];
	N = size_info[1];
	max_kd = (N > 200 ? N : 200);
	*one_delta = *two_delta = *trace_hat_out = 0;

	x = (double *) malloc(D * N * sizeof(double));
	x_tmp = (double *) malloc(D * N * sizeof(double));
	temp = (double *) malloc(N * sizeof(double));
	a_tmp = (long *) malloc(max_kd * sizeof(long));
	xi_tmp = (double *) malloc(max_kd * sizeof(double));
	vert_tmp = (double *) malloc(D * 2 * sizeof(double));
	vval_tmp = (double *) malloc((D + 1) * max_kd * sizeof(double));
	diag_tmp = (double *) malloc(N * sizeof(double));
	param_tmp = (long *) malloc(N * sizeof(long));
	order_parametric = (long *) malloc(D * sizeof(long));
	order_drop_sqr = (long *) malloc(D * sizeof(long));
        if((*iterations) > 0)
                pseudo_resid = (double *) malloc(N * sizeof(double));

	long* ltemp = (long*) temp;

	new_cell = (*span) * (*cell);
	for(i = 0; i < N; i++) 
		robust[i] = 1;
        for(i = 0; i < (N * D); i++)
                x_tmp[i] = x_[i];
	if((*normalize) && (D > 1)) {
		cut = (int)ceil(0.100000000000000000001 * N);
		for(i = 0; i < D; i++) {
			k = i * N;
			for(j = 0; j < N; j++)
				temp[j] = x_[k + j];
			qsort(temp, N, sizeof(double), vcomp);
			sum = 0;
			for(j = cut; j <= (N - cut - 1); j++)
			        sum = sum + temp[j];
			mean = sum / (N - 2 * cut);
			sum = 0;
			for(j = cut; j <= (N - cut - 1); j++) {
				temp[j] = temp[j] - mean;
				sum = sum + temp[j] * temp[j];
			}
			divisor[i] = sqrt(sum / (N - 2 * cut - 1));
			for(j = 0; j < N; j++) {
				p = k + j;
				x_tmp[p] = x_[p] / divisor[i];		
			}
		}
	}
	else
		for(i = 0; i < D; i++) divisor[i] = 1;
	j = D - 1;
	for(i = 0; i < D; i++) {
		sum_drop_sqr = sum_drop_sqr + drop_square[i];
		sum_parametric = sum_parametric + parametric[i];
		if(parametric[i])
			order_parametric[j--] = i;
		else
			order_parametric[nonparametric++] = i;
	}
        for(i = 0; i < D; i++) {
                order_drop_sqr[i] = 2 - drop_square[order_parametric[i]];
		k = i * N;
		p = order_parametric[i] * N;
	        for(j = 0; j < N; j++)
		        x[k + j] = x_tmp[p + j];
        }
	if((*degree) == 1 && sum_drop_sqr) {
		fprintf(stderr, "Specified the square of a factor predictor to be dropped when degree = 1");
		exit(1);
	}
	if(D == 1 && sum_drop_sqr) {
		fprintf(stderr, "Specified the square of a predictor to be dropped with only one numeric predictor");
		exit(1);
	}
	if(sum_parametric == D) {
		fprintf(stderr, "Specified parametric for all predictors");
		exit(1);
        }
	for(j = 0; j <= (*iterations); j++) {
		new_stat = j ? "none" : *statistics;
		for(i = 0; i < N; i++)
			robust[i] = weights[i] * robust[i];
		condition(surface, new_stat, trace_hat_in);
		setLf = !strcmp(surf_stat, "interpolate/exact");
		loess_raw(y, x, weights, robust, &D, &N, span, degree, 
			&nonparametric, order_drop_sqr, &sum_drop_sqr, 
			&new_cell, &surf_stat, fitted_values, parameter, a, 
			xi, vert, vval, diagonal, &trL, &delta1, &delta2, 
			&setLf); 
		if(j == 0) {
			*trace_hat_out = trL;
			*one_delta = delta1;
			*two_delta = delta2;
		}
		for(i = 0; i < N; i++)
			fitted_residuals[i] = y[i] - fitted_values[i];
		if(j < (*iterations))
			lowesw(fitted_residuals, &N, robust, ltemp);
	}
	if((*iterations) > 0) {
		lowesp(&N, y, fitted_values, weights, robust, ltemp, pseudovalues);
		
		loess_raw(pseudovalues, x, weights, weights, &D, &N, span, 
			degree,	&nonparametric, order_drop_sqr, &sum_drop_sqr,
			&new_cell, &surf_stat, temp, param_tmp, a_tmp, xi_tmp,
			vert_tmp, vval_tmp, diag_tmp, &trL_tmp, &d1_tmp, &d2_tmp, &zero);
		for(i = 0; i < N; i++)
			pseudo_resid[i] = pseudovalues[i] - temp[i];
	}
	if((*iterations) == 0)
		for(i = 0; i < N; i++)
			sum_squares = sum_squares + weights[i] * 
					fitted_residuals[i] * fitted_residuals[i];
	else 
		for(i = 0; i < N; i++)
			sum_squares = sum_squares + weights[i] *
					pseudo_resid[i] * pseudo_resid[i];
	*enp = (*one_delta) + 2 * (*trace_hat_out) - N;
	*s = sqrt(sum_squares / (*one_delta));

	free(x);
	free(x_tmp);
	free(temp);
	free(xi_tmp);
	free(vert_tmp);
	free(vval_tmp);
	free(diag_tmp);
	free(a_tmp);
	free(param_tmp);
	free(order_parametric);
	free(order_drop_sqr);
        if((*iterations) > 0)
                free(pseudo_resid);
}

void loess(loess_struct *lo)
{
	long	size_info[2], iterations;

	size_info[0] = lo->in.p;
	size_info[1] = lo->in.n;

	iterations = (!strcmp(lo->model.family, "gaussian")) ? 0 :
		lo->control.iterations;		
        if(!strcmp(lo->control.trace_hat, "wait.to.decide")) {
                if(!strcmp(lo->control.surface, "interpolate"))
                        lo->control.trace_hat = (lo->in.n < 500) ? "exact" : "approximate";
	        else 
		        lo->control.trace_hat = "exact";
        }
	loess_(lo->in.y, lo->in.x, size_info, lo->in.weights, 
		&lo->model.span,
		&lo->model.degree,
		lo->model.parametric,
		lo->model.drop_square,
		&lo->model.normalize,
		&lo->control.statistics,
		&lo->control.surface,
		&lo->control.cell,
		&lo->control.trace_hat,
		&iterations,
		lo->out.fitted_values,
		lo->out.fitted_residuals,
		&lo->out.enp,
		&lo->out.s,
		&lo->out.one_delta,
		&lo->out.two_delta,
		lo->out.pseudovalues,
		&lo->out.trace_hat,
		lo->out.diagonal,
		lo->out.robust,
		lo->out.divisor,
		lo->kd_tree.parameter,
		lo->kd_tree.a,
		lo->kd_tree.xi,
		lo->kd_tree.vert,
		lo->kd_tree.vval);
}

void loess_free_mem(struct	loess_struct *lo)
{
    free(lo->in.x);
	free(lo->in.y);
	free(lo->in.weights);
	free(lo->out.fitted_values);
	free(lo->out.fitted_residuals);
	free(lo->out.pseudovalues);
	free(lo->out.diagonal);
	free(lo->out.robust);
	free(lo->out.divisor);
	free(lo->kd_tree.parameter);
	free(lo->kd_tree.a);
	free(lo->kd_tree.xi);
	free(lo->kd_tree.vert);
	free(lo->kd_tree.vval);
}

void loess_summary(struct	loess_struct	*lo)
{
        printf("Number of Observations: %ld\n", lo->in.n);
	printf("Equivalent Number of Parameters: %.1f\n", lo->out.enp);
	if(!strcmp(lo->model.family, "gaussian"))
		printf("Residual Standard Error: ");
	else
		printf("Residual Scale Estimate: ");
	printf("%.4f\n", lo->out.s);
}

// Loess Internal functions

static long c__2 = 2;
static long c__180 = 180;
static long c__4 = 4;
static long c__1 = 1;
static long c__120 = 120;
static long c__121 = 121;
static long c__0 = 0;
static long c__1000 = 1000;
static long c__15 = 15;
static long c__21 = 21;
static long c__182 = 182;
static long c__101 = 101;
static long c__193 = 193;
static long c__181 = 181;
static long c__122 = 122;
static long c__104 = 104;
static long c__105 = 105;
static long c__123 = 123;
static long c__10000 = 10000;
static long c__194 = 194;
static long c__196 = 196;
static long c__174 = 174;
static long c__171 = 171;
static long c__100 = 100;
static long c__195 = 195;
static long c__102 = 102;
static long c__103 = 103;
static long c__172 = 172;
static long c__173 = 173;
static long c__186 = 186;
static long c__175 = 175;
static long c__187 = 187;
static long c__185 = 185;

/* Subroutine */ int ehg126_(long *d, long *n, long *vc, double *
	x, double *v, long *nvmax)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long v_dim1, v_offset, x_dim1, x_offset, i__1, i__2;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    static double beta;
    static long i, j, k;
    static double t, alpha;
    extern double d1mach_(long *);
    static double machin, mu;

    /* Parameter adjustments */
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
/*     MachInf -> machin */
    ++execnt;
    if (execnt == 1) {
	machin = d1mach_(&c__2);
    }
/*     fill in vertices for bounding box of $x$ */
/*     lower left, upper right */
    i__1 = *d;
    for (k = 1; k <= i__1; ++k) {
	alpha = machin;
	beta = -machin;
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    t = x[i + k * x_dim1];
	    alpha = min(alpha,t);
	    beta = max(beta,t);
/* L4: */
	}
/*        expand the box a little */
/* Computing MAX */
/* Computing MAX */
	d__3 = abs(alpha), d__4 = abs(beta);
	d__1 = beta - alpha, d__2 = max(d__3,d__4) * 1e-10 + 1e-30;
	mu = max(d__1,d__2) * .005;
	alpha -= mu;
	beta += mu;
	v[k * v_dim1 + 1] = alpha;
	v[*vc + k * v_dim1] = beta;
/* L3: */
    }
/*     remaining vertices */
    i__1 = *vc - 1;
    for (i = 2; i <= i__1; ++i) {
	j = i - 1;
	i__2 = *d;
	for (k = 1; k <= i__2; ++k) {
	    v[i + k * v_dim1] = v[j % 2 * (*vc - 1) + 1 + k * v_dim1];
	    j = (long) ((double) j / 2.);
/* L6: */
	}
/* L5: */
    }
    return 0;
} /* ehg126_ */

/* Subroutine */ int ehg125_(long *p, long *nv, double *v, long *
	vhit, long *nvmax, long *d, long *k, double *t, long *
	r, long *s, long *f, long *l, long *u)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long f_dim1, f_offset, l_dim1, l_offset, u_dim1, u_offset, v_dim1, 
	    v_offset, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int ehg182_(long *);
    static long h, i, j, m;
    static long match, i1, i2;
    static long i3, mm;

    /* Parameter adjustments */
    u_dim1 = *r;
    u_offset = (u_dim1 << 1) + 1;
    u -= u_offset;
    l_dim1 = *r;
    l_offset = (l_dim1 << 1) + 1;
    l -= l_offset;
    f_dim1 = *r;
    f_offset = (f_dim1 << 1) + 1;
    f -= f_offset;
    --vhit;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;

    /* Function Body */
    ++execnt;
    h = *nv;
    i__1 = *r;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *s;
	for (j = 1; j <= i__2; ++j) {
	    ++h;
	    i__3 = *d;
	    for (i3 = 1; i3 <= i__3; ++i3) {
		v[h + i3 * v_dim1] = v[f[i + (j << 1) * f_dim1] + i3 * v_dim1]
			;
/* L5: */
	    }
	    v[h + *k * v_dim1] = *t;
/*           check for redundant vertex */
	    match = false;
	    m = 1;
/*           top of while loop */
L6:
	    if (! match) {
		i1 = m <= *nv;
	    } else {
		i1 = false;
	    }
	    if (! i1) {
		goto L7;
	    }
	    match = v[m + v_dim1] == v[h + v_dim1];
	    mm = 2;
/*              top of while loop */
L8:
	    if (match) {
		i2 = mm <= *d;
	    } else {
		i2 = false;
	    }
	    if (! i2) {
		goto L9;
	    }
	    match = v[m + mm * v_dim1] == v[h + mm * v_dim1];
	    ++mm;
	    goto L8;
/*              bottom of while loop */
L9:
	    ++m;
	    goto L6;
/*           bottom of while loop */
L7:
	    --m;
	    if (match) {
		--h;
	    } else {
		m = h;
		if (vhit[1] >= 0) {
		    vhit[m] = *p;
		}
	    }
	    l[i + (j << 1) * l_dim1] = f[i + (j << 1) * f_dim1];
	    l[i + ((j << 1) + 1) * l_dim1] = m;
	    u[i + (j << 1) * u_dim1] = m;
	    u[i + ((j << 1) + 1) * u_dim1] = f[i + ((j << 1) + 1) * f_dim1];
/* L4: */
	}
/* L3: */
    }
    *nv = h;
    if (! (*nv <= *nvmax)) {
	ehg182_(&c__180);
    }
    return 0;
} /* ehg125_ */

long ehg138_(long *i, double *z, long *a, double *xi, 
	long *lo, long *hi, long *ncmax)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long ret_val;

    /* Local variables */
    static long j;
    static long i1;

    /* Parameter adjustments */
    --hi;
    --lo;
    --xi;
    --a;
    --z;

    /* Function Body */
    ++execnt;
/*     descend tree until leaf or ambiguous */
    j = *i;
/*     top of while loop */
L3:
    if (a[j] != 0) {
	i1 = z[a[j]] != xi[j];
    } else {
	i1 = false;
    }
    if (! i1) {
	goto L4;
    }
    if (z[a[j]] < xi[j]) {
	j = lo[j];
    } else {
	j = hi[j];
    }
    goto L3;
/*     bottom of while loop */
L4:
    ret_val = j;
    return ret_val;
} /* ehg138_ */

/* Subroutine */ int ehg106_(long *il, long *ir, long *k, long *
	nk, double *p, long *pi, long *n)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long p_dim1, p_offset;

    /* Local variables */
    static long i, j, l, r;
    static double t;
    static long ii;

    /* Parameter adjustments */
    --pi;
    p_dim1 = *nk;
    p_offset = p_dim1 + 1;
    p -= p_offset;

    /* Function Body */
    ++execnt;
/*     find the $k$-th smallest of $n$ elements */
/*     Floyd+Rivest, CACM Mar '75, Algorithm 489 */
    l = *il;
    r = *ir;
/*     top of while loop */
L3:
    if (! (l < r)) {
	goto L4;
    }
/*        to avoid recursion, sophisticated partition deleted */
/*        partition $x sub {l..r}$ about $t$ */
    t = p[pi[*k] * p_dim1 + 1];
    i = l;
    j = r;
    ii = pi[l];
    pi[l] = pi[*k];
    pi[*k] = ii;
    if (t < p[pi[r] * p_dim1 + 1]) {
	ii = pi[l];
	pi[l] = pi[r];
	pi[r] = ii;
    }
/*        top of while loop */
L5:
    if (! (i < j)) {
	goto L6;
    }
    ii = pi[i];
    pi[i] = pi[j];
    pi[j] = ii;
    ++i;
    --j;
/*           top of while loop */
L7:
    if (! (p[pi[i] * p_dim1 + 1] < t)) {
	goto L8;
    }
    ++i;
    goto L7;
/*           bottom of while loop */
L8:
/*           top of while loop */
L9:
    if (! (t < p[pi[j] * p_dim1 + 1])) {
	goto L10;
    }
    --j;
    goto L9;
/*           bottom of while loop */
L10:
    goto L5;
/*        bottom of while loop */
L6:
    if (p[pi[l] * p_dim1 + 1] == t) {
	ii = pi[l];
	pi[l] = pi[j];
	pi[j] = ii;
    } else {
	++j;
	ii = pi[r];
	pi[r] = pi[j];
	pi[j] = ii;
    }
    if (j <= *k) {
	l = j + 1;
    }
    if (*k <= j) {
	r = j - 1;
    }
    goto L3;
/*     bottom of while loop */
L4:
    return 0;
} /* ehg106_ */

/* Subroutine */ int ehg127_(double *q, long *n, long *d, long *
	nf, double *f, double *x, long *psi, double *y, 
	double *rw, long *kernel, long *k, double *dist, 
	double *eta, double *b, long *od, double *w, 
	double *rcond, long *sing, double *sigma, double *u, 
	double *e, double *dgamma, double *qraux, double *
	work, double *tol, long *dd, long *tdeg, long *cdeg, 
	double *s)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static double scal;
    extern double ddot_(long *, double *, long *, double *, 
	    long *);
    static long info, jpvt;
    extern /* Subroutine */ int ehg106_(long *, long *, long *, 
	    long *, double *, long *, long *), ehg182_(long *)
	    , ehg184_(char *, double *, long *, long *, long);
    static double g[15];
    static long i, j;
    extern /* Subroutine */ int dqrdc_(double *, long *, long *, 
	    long *, double *, long *, double *, long *), 
	    dsvdc_(double *, long *, long *, long *, double *
	    , double *, double *, long *, double *, long *, 
	    double *, long *, long *), dqrsl_(double *, long 
	    *, long *, long *, double *, double *, double *,
	     double *, double *, double *, double *, long *
	    , long *);
    static double i1, i2;
    static long i3;
    static double i4, i5, i6, i7;
    static long i9;
    static double i8;
    extern double d1mach_(long *);
    static long inorm2;
    static double i10;
    static long jj;
    static double machep;
    extern long idamax_(long *, double *, long *);
    static double colnor[15];
    static long column;
    static double rho;

    /* Parameter adjustments */
    --cdeg;
    --work;
    --qraux;
    --dgamma;
    e -= 16;
    u -= 16;
    --sigma;
    --w;
    b_dim1 = *nf;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    --eta;
    --dist;
    --rw;
    --y;
    --psi;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --q;

    /* Function Body */
/*     colnorm -> colnor */
/*     E -> g */
/*     MachEps -> machep */
/*     V -> e */
/*     X -> b */
    ++execnt;
    if (execnt == 1) {
	machep = d1mach_(&c__4);
    }
/*     sort by distance */
    i__1 = *n;
    for (i3 = 1; i3 <= i__1; ++i3) {
	dist[i3] = 0.;
/* L3: */
    }
    i__1 = *dd;
    for (j = 1; j <= i__1; ++j) {
	i4 = q[j];
	i__2 = *n;
	for (i3 = 1; i3 <= i__2; ++i3) {
/* Computing 2nd power */
	    d__1 = x[i3 + j * x_dim1] - i4;
	    dist[i3] += d__1 * d__1;
/* L5: */
	}
/* L4: */
    }
    ehg106_(&c__1, n, nf, &c__1, &dist[1], &psi[1], n);
    rho = dist[psi[*nf]] * max(1.,*f);
    if (! (0. < rho)) {
	ehg182_(&c__120);
    }
/*     compute neighborhood weights */
    if (*kernel == 2) {
	i__1 = *nf;
	for (i = 1; i <= i__1; ++i) {
	    if (dist[psi[i]] < rho) {
		i1 = sqrt(rw[psi[i]]);
	    } else {
		i1 = 0.;
	    }
	    w[i] = i1;
/* L6: */
	}
    } else {
	i__1 = *nf;
	for (i3 = 1; i3 <= i__1; ++i3) {
	    w[i3] = sqrt(dist[psi[i3]] / rho);
/* L7: */
	}
	i__1 = *nf;
	for (i3 = 1; i3 <= i__1; ++i3) {
/* Computing 3rd power */
	    d__2 = w[i3], d__3 = d__2;
/* Computing 3rd power */
	    d__1 = 1 - d__3 * (d__2 * d__2), d__4 = d__1;
	    w[i3] = sqrt(rw[psi[i3]] * (d__4 * (d__1 * d__1)));
/* L8: */
	}
    }
    if ((d__1 = w[idamax_(nf, &w[1], &c__1)], abs(d__1)) == 0.) {
	ehg184_("at ", &q[1], dd, &c__1, 3L);
	ehg184_("radius ", &rho, &c__1, &c__1, 7L);
	if (true) {
	    ehg182_(&c__121);
	}
    }
/*     fill design matrix */
    column = 1;
    i__1 = *nf;
    for (i3 = 1; i3 <= i__1; ++i3) {
	b[i3 + column * b_dim1] = w[i3];
/* L9: */
    }
    if (*tdeg >= 1) {
	i__1 = *d;
	for (j = 1; j <= i__1; ++j) {
	    if (cdeg[j] >= 1) {
		++column;
		i5 = q[j];
		i__2 = *nf;
		for (i3 = 1; i3 <= i__2; ++i3) {
		    b[i3 + column * b_dim1] = w[i3] * (x[psi[i3] + j * x_dim1]
			     - i5);
/* L11: */
		}
	    }
/* L10: */
	}
    }
    if (*tdeg >= 2) {
	i__1 = *d;
	for (j = 1; j <= i__1; ++j) {
	    if (cdeg[j] >= 1) {
		if (cdeg[j] >= 2) {
		    ++column;
		    i6 = q[j];
		    i__2 = *nf;
		    for (i3 = 1; i3 <= i__2; ++i3) {
/* Computing 2nd power */
			d__1 = x[psi[i3] + j * x_dim1] - i6;
			b[i3 + column * b_dim1] = w[i3] * (d__1 * d__1);
/* L13: */
		    }
		}
		i__2 = *d;
		for (jj = j + 1; jj <= i__2; ++jj) {
		    if (cdeg[jj] >= 1) {
			++column;
			i7 = q[j];
			i8 = q[jj];
			i__3 = *nf;
			for (i3 = 1; i3 <= i__3; ++i3) {
			    b[i3 + column * b_dim1] = w[i3] * (x[psi[i3] + j *
				     x_dim1] - i7) * (x[psi[i3] + jj * x_dim1]
				     - i8);
/* L15: */
			}
		    }
/* L14: */
		}
	    }
/* L12: */
	}
	*k = column;
    }
    i__1 = *nf;
    for (i3 = 1; i3 <= i__1; ++i3) {
	eta[i3] = w[i3] * y[psi[i3]];
/* L16: */
    }
/*     equilibrate columns */
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	scal = 0.;
	i__2 = *nf;
	for (inorm2 = 1; inorm2 <= i__2; ++inorm2) {
/* Computing 2nd power */
	    d__1 = b[inorm2 + j * b_dim1];
	    scal += d__1 * d__1;
/* L18: */
	}
	scal = sqrt(scal);
	if (0. < scal) {
	    i__2 = *nf;
	    for (i3 = 1; i3 <= i__2; ++i3) {
		b[i3 + j * b_dim1] /= scal;
/* L19: */
	    }
	    colnor[j - 1] = scal;
	} else {
	    colnor[j - 1] = 1.;
	}
/* L17: */
    }
/*     singular value decomposition */
    dqrdc_(&b[b_offset], nf, nf, k, &qraux[1], &jpvt, &work[1], &c__0);
    dqrsl_(&b[b_offset], nf, nf, k, &qraux[1], &eta[1], &work[1], &eta[1], &
	    eta[1], &work[1], &work[1], &c__1000, &info);
    i__1 = *k;
    for (i9 = 1; i9 <= i__1; ++i9) {
	i__2 = *k;
	for (i3 = 1; i3 <= i__2; ++i3) {
	    u[i3 + i9 * 15] = 0.;
/* L21: */
	}
/* L20: */
    }
    i__1 = *k;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *k;
	for (j = i; j <= i__2; ++j) {
	    u[i + j * 15] = b[i + j * b_dim1];
/* L23: */
	}
/* L22: */
    }
    dsvdc_(&u[16], &c__15, k, k, &sigma[1], g, &u[16], &c__15, &e[16], &c__15,
	     &work[1], &c__21, &info);
    if (! (info == 0)) {
	ehg182_(&c__182);
    }
    *tol = sigma[1] * (machep * 100);
/* Computing MIN */
    d__1 = *rcond, d__2 = sigma[*k] / sigma[1];
    *rcond = min(d__1,d__2);
    if (sigma[*k] <= *tol) {
	++(*sing);
	if (*sing == 1) {
	    ehg184_("Warning. pseudoinverse used at", &q[1], d, &c__1, 30L);
	    d__1 = sqrt(rho);
	    ehg184_("neighborhood radius", &d__1, &c__1, &c__1, 19L);
	    ehg184_("reciprocal condition number ", rcond, &c__1, &c__1, 28L);
	} else {
	    if (*sing == 2) {
		ehg184_("There are other near singularities as well.", &rho, &
			c__1, &c__1, 43L);
	    }
	}
    }
/*     compensate for equilibration */
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i10 = colnor[j - 1];
	i__2 = *k;
	for (i3 = 1; i3 <= i__2; ++i3) {
	    e[j + i3 * 15] /= i10;
/* L25: */
	}
/* L24: */
    }
/*     solve least squares problem */
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	if (*tol < sigma[j]) {
	    i2 = ddot_(k, &u[j * 15 + 1], &c__1, &eta[1], &c__1) / sigma[j];
	} else {
	    i2 = 0.;
	}
	dgamma[j] = i2;
/* L26: */
    }
    i__1 = *od;
    for (j = 0; j <= i__1; ++j) {
	s[j] = ddot_(k, &e[j + 16], &c__15, &dgamma[1], &c__1);
/* L27: */
    }
    return 0;
} /* ehg127_ */

/* Subroutine */ int ehg131_(double *x, double *y, double *rw, 
	double *trl, double *diagl, long *kernel, long *k, 
	long *n, long *d, long *nc, long *ncmax, long *vc, 
	long *nv, long *nvmax, long *nf, double *f, long *a, 
	long *c, long *hi, long *lo, long *pi, long *psi, 
	double *v, long *vhit, double *vval, double *xi, 
	double *dist, double *eta, double *b, long *ntol, 
	double *fd, double *w, double *vval2, double *rcond, 
	long *sing, long *dd, long *tdeg, long *cdeg, long *lq,
	 double *lf, long *setlf)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long lq_dim1, lq_offset, c_dim1, c_offset, lf_dim1, lf_dim2, lf_offset,
	     v_dim1, v_offset, vval_dim1, vval_offset, vval2_dim1, 
	    vval2_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int ehg124_(long *, long *, long *, 
	    long *, long *, long *, long *, long *, double 
	    *, long *, long *, double *, long *, long *, 
	    long *, double *, long *, long *, long *, 
	    double *, long *), ehg126_(long *, long *, long *,
	     double *, double *, long *), ehg182_(long *), 
	    ehg139_(double *, long *, long *, long *, long *, 
	    long *, double *, double *, long *, long *, 
	    double *, double *, double *, long *, long *, 
	    double *, long *, double *, double *, long *,
	     double *, double *, double *, long *, long *, 
	    long *, double *, long *, long *, long *, long 
	    *, double *, long *, long *, long *, long *, 
	    long *, double *, long *, double *);
    extern double dnrm2_(long *, double *, long *);
    static long j;
    static double delta[8];
    static long i1, i2, identi;

    /* Parameter adjustments */
    lf_dim1 = *d + 1;
    lf_dim2 = *nvmax;
    lf_offset = lf_dim1 * (lf_dim2 + 1);
    lf -= lf_offset;
    lq_dim1 = *nvmax;
    lq_offset = lq_dim1 + 1;
    lq -= lq_offset;
    --cdeg;
    vval2_dim1 = *d + 1;
    vval2_offset = vval2_dim1;
    vval2 -= vval2_offset;
    --w;
    --b;
    --eta;
    --dist;
    --xi;
    vval_dim1 = *d + 1;
    vval_offset = vval_dim1;
    vval -= vval_offset;
    --vhit;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --psi;
    --pi;
    --lo;
    --hi;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --a;
    --diagl;
    --rw;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
/*     Identity -> identi */
/*     X -> b */
    ++execnt;
    if (! (*d <= 8)) {
	ehg182_(&c__101);
    }
/*     build $k$-d tree */
    ehg126_(d, n, vc, &x[x_offset], &v[v_offset], nvmax);
    *nv = *vc;
    *nc = 1;
    i__1 = *vc;
    for (j = 1; j <= i__1; ++j) {
	c[j + *nc * c_dim1] = j;
	vhit[j] = 0;
/* L3: */
    }
    i__1 = *d;
    for (i1 = 1; i1 <= i__1; ++i1) {
	delta[i1 - 1] = v[*vc + i1 * v_dim1] - v[i1 * v_dim1 + 1];
/* L4: */
    }
    *fd *= dnrm2_(d, delta, &c__1);
    i__1 = *n;
    for (identi = 1; identi <= i__1; ++identi) {
	pi[identi] = identi;
/* L5: */
    }
    ehg124_(&c__1, n, d, n, nv, nc, ncmax, vc, &x[x_offset], &pi[1], &a[1], &
	    xi[1], &lo[1], &hi[1], &c[c_offset], &v[v_offset], &vhit[1], 
	    nvmax, ntol, fd, dd);
/*     smooth */
    if (*trl != 0.) {
	i__1 = *nv;
	for (i2 = 1; i2 <= i__1; ++i2) {
	    i__2 = *d;
	    for (i1 = 0; i1 <= i__2; ++i1) {
		vval2[i1 + i2 * vval2_dim1] = 0.;
/* L7: */
	    }
/* L6: */
	}
    }
    ehg139_(&v[v_offset], nvmax, nv, n, d, nf, f, &x[x_offset], &pi[1], &psi[
	    1], &y[1], &rw[1], trl, kernel, k, &dist[1], (long*)&dist[1], &eta[1], &
	    b[1], d, &w[1], &diagl[1], &vval2[vval2_offset], nc, vc, &a[1], &
	    xi[1], &lo[1], &hi[1], &c[c_offset], &vhit[1], rcond, sing, dd, 
	    tdeg, &cdeg[1], &lq[lq_offset], &lf[lf_offset], setlf, &vval[
	    vval_offset]);
    return 0;
} /* ehg131_ */

/* Subroutine */ int ehg133_(long *n, long *d, long *vc, long *
	nvmax, long *nc, long *ncmax, long *a, long *c, long *
	hi, long *lo, double *v, double *vval, double *xi, 
	long *m, double *z, double *s)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long c_dim1, c_offset, v_dim1, v_offset, vval_dim1, vval_offset, 
	    z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    extern double ehg128_(double *, long *, long *, long *, 
	    long *, double *, long *, long *, long *, 
	    double *, long *, double *);
    static long i;
    static double delta[8];
    static long i1;

    /* Parameter adjustments */
    --s;
    z_dim1 = *m;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --xi;
    vval_dim1 = *d + 1;
    vval_offset = vval_dim1;
    vval -= vval_offset;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --lo;
    --hi;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --a;

    /* Function Body */
    ++execnt;
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *d;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    delta[i1 - 1] = z[i + i1 * z_dim1];
/* L4: */
	}
	s[i] = ehg128_(delta, d, ncmax, vc, &a[1], &xi[1], &lo[1], &hi[1], &c[
		c_offset], &v[v_offset], nvmax, &vval[vval_offset]);
/* L3: */
    }
    return 0;
} /* ehg133_ */

/* Subroutine */ int ehg140_(long *iw, long *i, long *j)
{
    /* Initialized data */

    static long execnt = 0;

    /* Parameter adjustments */
    --iw;

    /* Function Body */
    ++execnt;
    iw[*i] = *j;
    return 0;
} /* ehg140_ */

/* Subroutine */ int ehg141_(double *trl, long *n, long *deg, 
	long *k, long *d, long *nsing, long *dk, double *
	delta1, double *delta2)
{
    /* Initialized data */

    static double c[48] = { .297162,.380266,.5886043,.4263766,.3346498,
	    .6271053,.5241198,.3484836,.6687687,.6338795,.4076457,.7207693,
	    .1611761,.3091323,.4401023,.2939609,.3580278,.5555741,.397239,
	    .4171278,.6293196,.4675173,.469907,.6674802,.2848308,.2254512,
	    .2914126,.5393624,.251723,.389897,.7603231,.2969113,.474013,
	    .9664956,.3629838,.5348889,.207567,.2822574,.2369957,.3911566,
	    .2981154,.3623232,.5508869,.3501989,.4371032,.7002667,.4291632,
	    .493037 };

    /* System generated locals */
    double d__1, d__2;

    /* Builtin functions */
    double sqrt(double), exp(double), pow_dd(double *, double 
	    *);

    /* Local variables */
    static double corx;
    extern /* Subroutine */ int ehg184_(char *, double *, long *, 
	    long *, long);
    extern double ehg176_(double *);
    static long i;
    static double z, c1, c2, c3, c4;

/*     coef, d, deg, del */
    if (*deg == 0) {
	*dk = 1;
    }
    if (*deg == 1) {
	*dk = *d + 1;
    }
    if (*deg == 2) {
	*dk = (long) ((double) ((*d + 2) * (*d + 1)) / 2.);
    }
    corx = sqrt(*k / (double) (*n));
    z = (sqrt(*k / *trl) - corx) / (1 - corx);
    if (*nsing == 0 && 1. < z) {
	ehg184_("Chernobyl! trL<k", trl, &c__1, &c__1, 16L);
    }
    if (z < 0.) {
	ehg184_("Chernobyl! trL>n", trl, &c__1, &c__1, 16L);
    }
/* Computing MIN */
    d__1 = 1., d__2 = max(0.,z);
    z = min(d__1,d__2);
    c4 = exp(ehg176_(&z));
    i = (min(*d,4) - 1 + ((long)(*deg - 1) << 2)) * 3 + 1;
    if (*d <= 4) {
	c1 = c[i - 1];
	c2 = c[i];
	c3 = c[i + 1];
    } else {
	c1 = c[i - 1] + (*d - 4) * (c[i - 1] - c[i - 4]);
	c2 = c[i] + (*d - 4) * (c[i] - c[i - 3]);
	c3 = c[i + 1] + (*d - 4) * (c[i + 1] - c[i - 2]);
    }
    d__1 = 1 - z;
    *delta1 = *n - *trl * exp(c1 * pow_dd(&z, &c2) * pow_dd(&d__1, &c3) * c4);
    i += 24;
    if (*d <= 4) {
	c1 = c[i - 1];
	c2 = c[i];
	c3 = c[i + 1];
    } else {
	c1 = c[i - 1] + (*d - 4) * (c[i - 1] - c[i - 4]);
	c2 = c[i] + (*d - 4) * (c[i] - c[i - 3]);
	c3 = c[i + 1] + (*d - 4) * (c[i + 1] - c[i - 2]);
    }
    d__1 = 1 - z;
    *delta2 = *n - *trl * exp(c1 * pow_dd(&z, &c2) * pow_dd(&d__1, &c3) * c4);
    return 0;
} /* ehg141_ */

/* Subroutine */ int lowesc(long *n, double *l, double *ll, 
	double *trl, double *delta1, double *delta2)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long l_dim1, l_offset, ll_dim1, ll_offset, i__1, i__2;

    /* Local variables */
    extern double ddot_(long *, double *, long *, double *, 
	    long *);
    static long i, j;

    /* Parameter adjustments */
    ll_dim1 = *n;
    ll_offset = ll_dim1 + 1;
    ll -= ll_offset;
    l_dim1 = *n;
    l_offset = l_dim1 + 1;
    l -= l_offset;

    /* Function Body */
    ++execnt;
/*     compute $LL~=~(I-L)(I-L)'$ */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	--l[i + i * l_dim1];
/* L3: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = i;
	for (j = 1; j <= i__2; ++j) {
	    ll[i + j * ll_dim1] = ddot_(n, &l[i + l_dim1], n, &l[j + l_dim1], 
		    n);
/* L5: */
	}
/* L4: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *n;
	for (j = i + 1; j <= i__2; ++j) {
	    ll[i + j * ll_dim1] = ll[j + i * ll_dim1];
/* L7: */
	}
/* L6: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	++l[i + i * l_dim1];
/* L8: */
    }
/*     accumulate first two traces */
    *trl = 0.;
    *delta1 = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	*trl += l[i + i * l_dim1];
	*delta1 += ll[i + i * ll_dim1];
/* L9: */
    }
/*     $delta sub 2 = "tr" LL sup 2$ */
    *delta2 = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	*delta2 += ddot_(n, &ll[i + ll_dim1], n, &ll[i * ll_dim1 + 1], &c__1);
/* L10: */
    }
    return 0;
} /* lowesc_ */

/* Subroutine */ int ehg169_(long *d, long *vc, long *nc, long *
	ncmax, long *nv, long *nvmax, double *v, long *a, 
	double *xi, long *c, long *hi, long *lo)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long c_dim1, c_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    double d__1;

    /* Builtin functions */
    long pow_ii(long *, long *);

    /* Local variables */
    extern /* Subroutine */ int ehg125_(long *, long *, double *, 
	    long *, long *, long *, long *, double *, long 
	    *, long *, long *, long *, long *), ehg182_(long *)
	    ;
    static long i, j, k, p, mc, mv;
    extern long ifloor_(double *);
    static long novhit[1];

    /* Parameter adjustments */
    --lo;
    --hi;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --xi;
    --a;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;

    /* Function Body */
    ++execnt;
/*     as in bbox */
/*     remaining vertices */
    i__1 = *vc - 1;
    for (i = 2; i <= i__1; ++i) {
	j = i - 1;
	i__2 = *d;
	for (k = 1; k <= i__2; ++k) {
	    v[i + k * v_dim1] = v[j % 2 * (*vc - 1) + 1 + k * v_dim1];
	    d__1 = (double) j / 2.;
	    j = ifloor_(&d__1);
/* L4: */
	}
/* L3: */
    }
/*     as in ehg131 */
    mc = 1;
    mv = *vc;
    novhit[0] = -1;
    i__1 = *vc;
    for (j = 1; j <= i__1; ++j) {
	c[j + mc * c_dim1] = j;
/* L5: */
    }
/*     as in rbuild */
    p = 1;
/*     top of while loop */
L6:
    if (! (p <= *nc)) {
	goto L7;
    }
    if (a[p] != 0) {
	k = a[p];
/*           left son */
	++mc;
	lo[p] = mc;
/*           right son */
	++mc;
	hi[p] = mc;
	i__2 = k - 1;
	i__1 = pow_ii(&c__2, &i__2);
	i__4 = *d - k;
	i__3 = pow_ii(&c__2, &i__4);
	ehg125_(&p, &mv, &v[v_offset], novhit, nvmax, d, &k, &xi[p], &i__1, &
		i__3, &c[p * c_dim1 + 1], &c[lo[p] * c_dim1 + 1], &c[hi[p] * 
		c_dim1 + 1]);
    }
    ++p;
    goto L6;
/*     bottom of while loop */
L7:
    if (! (mc == *nc)) {
	ehg182_(&c__193);
    }
    if (! (mv == *nv)) {
	ehg182_(&c__193);
    }
    return 0;
} /* ehg169_ */

double ehg176_(double *z)
{
    /* Initialized data */

    static long d = 1;
    static long vc = 2;
    static long nv = 10;
    static long nc = 17;
    static long a[17] = { 1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0 };
    static struct {
	long e_1[7];
	long fill_2[7];
	long e_3;
	long fill_4[2];
	} equiv_94 = { 3, 5, 7, 9, 11, 13, 15, {0}, 17 };

#define hi ((long *)&equiv_94)

    static struct {
	long e_1[7];
	long fill_2[7];
	long e_3;
	long fill_4[2];
	} equiv_95 = { 2, 4, 6, 8, 10, 12, 14, {0}, 16 };

#define lo ((long *)&equiv_95)

    static struct {
	double e_1[7];
	double fill_2[7];
	double e_3;
	double fill_4[2];
	} equiv_96 = { .3705, .2017, .5591, .1204, .2815, .4536, .7132, {0}, 
		.8751 };

#define xi ((double *)&equiv_96)

    static long c[34]	/* was [2][17] */ = { 1,2,1,3,3,2,1,4,4,3,3,5,
	    5,2,1,6,6,4,4,7,7,3,3,8,8,5,5,9,9,2,9,10,10,2 };
    static double vval[20]	/* was [2][10] */ = { -.090572,4.4844,
	    -.010856,-.7736,-.053718,-.3495,.026152,-.7286,-.058387,.1611,
	    .095807,-.7978,-.031926,-.4457,-.06417,.032813,-.020636,.335,
	    .040172,-.041032 };
    static double v[10]	/* was [10][1] */ = { -.005,1.005,.3705,.2017,
	    .5591,.1204,.2815,.4536,.7132,.8751 };

    /* System generated locals */
    double ret_val;

    /* Local variables */
    extern double ehg128_(double *, long *, long *, long *, 
	    long *, double *, long *, long *, long *, 
	    double *, long *, double *);

    /* Parameter adjustments */
    --z;

    /* Function Body */
    ret_val = ehg128_(&z[1], &d, &nc, &vc, a, xi, lo, hi, c, v, &nv, vval);
    return ret_val;
} /* ehg176_ */

#undef xi
#undef lo
#undef hi


/* Subroutine */ int lowesa(double *trl, long *n, long *d, long 
	*tau, long *nsing, double *delta1, double *delta2)
{
    /* Initialized data */

    static long execnt = 0;

    extern /* Subroutine */ int ehg141_(double *, long *, long *, 
	    long *, long *, long *, long *, double *, 
	    double *);
    static double alpha, d1a, d1b, d2a, d2b;
    static long dka, dkb;

    ++execnt;
    ehg141_(trl, n, &c__1, tau, d, nsing, &dka, &d1a, &d2a);
    ehg141_(trl, n, &c__2, tau, d, nsing, &dkb, &d1b, &d2b);
    alpha = (double) (*tau - dka) / (double) (dkb - dka);
    *delta1 = (1 - alpha) * d1a + alpha * d1b;
    *delta2 = (1 - alpha) * d2a + alpha * d2b;
    return 0;
} /* lowesa_ */

/* Subroutine */ int ehg191_(long *m, double *z, double *l, 
	long *d, long *n, long *nf, long *nv, long *ncmax, 
	long *vc, long *a, double *xi, long *lo, long *hi, 
	long *c, double *v, long *nvmax, double *vval2, 
	double *lf, long *lq)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long lq_dim1, lq_offset, c_dim1, c_offset, l_dim1, l_offset, lf_dim1, 
	    lf_dim2, lf_offset, v_dim1, v_offset, vval2_dim1, vval2_offset, 
	    z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    extern double ehg128_(double *, long *, long *, long *, 
	    long *, double *, long *, long *, long *, 
	    double *, long *, double *);
    static long i, j, p, i1, i2;
    static double zi[8];
    static long lq1;

    /* Parameter adjustments */
    lq_dim1 = *nvmax;
    lq_offset = lq_dim1 + 1;
    lq -= lq_offset;
    lf_dim1 = *d + 1;
    lf_dim2 = *nvmax;
    lf_offset = lf_dim1 * (lf_dim2 + 1);
    lf -= lf_offset;
    vval2_dim1 = *d + 1;
    vval2_offset = vval2_dim1;
    vval2 -= vval2_offset;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --hi;
    --lo;
    --xi;
    --a;
    l_dim1 = *m;
    l_offset = l_dim1 + 1;
    l -= l_offset;
    z_dim1 = *m;
    z_offset = z_dim1 + 1;
    z -= z_offset;

    /* Function Body */
    ++execnt;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nv;
	for (i2 = 1; i2 <= i__2; ++i2) {
	    i__3 = *d;
	    for (i1 = 0; i1 <= i__3; ++i1) {
		vval2[i1 + i2 * vval2_dim1] = 0.;
/* L5: */
	    }
/* L4: */
	}
	i__2 = *nv;
	for (i = 1; i <= i__2; ++i) {
/*           linear search for i in Lq */
	    lq1 = lq[i + lq_dim1];
	    lq[i + lq_dim1] = j;
	    p = *nf;
/*           top of while loop */
L7:
	    if (! (lq[i + p * lq_dim1] != j)) {
		goto L8;
	    }
	    --p;
	    goto L7;
/*           bottom of while loop */
L8:
	    lq[i + lq_dim1] = lq1;
	    if (lq[i + p * lq_dim1] == j) {
		i__3 = *d;
		for (i1 = 0; i1 <= i__3; ++i1) {
		    vval2[i1 + i * vval2_dim1] = lf[i1 + (i + p * lf_dim2) * 
			    lf_dim1];
/* L9: */
		}
	    }
/* L6: */
	}
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
	    i__3 = *d;
	    for (i1 = 1; i1 <= i__3; ++i1) {
		zi[i1 - 1] = z[i + i1 * z_dim1];
/* L11: */
	    }
	    l[i + j * l_dim1] = ehg128_(zi, d, ncmax, vc, &a[1], &xi[1], &lo[
		    1], &hi[1], &c[c_offset], &v[v_offset], nvmax, &vval2[
		    vval2_offset]);
/* L10: */
	}
/* L3: */
    }
    return 0;
} /* ehg191_ */

/* Subroutine */ int ehg196(long *tau, long *d, double *f, 
	double *trl)
{
    /* Initialized data */

    static long execnt = 0;

    static double trla, trlb;
    extern /* Subroutine */ int ehg197_(long *, long *, long *, 
	    double *, long *, double *);
    static double alpha;
    static long dka, dkb;

    ++execnt;
    ehg197_(&c__1, tau, d, f, &dka, &trla);
    ehg197_(&c__2, tau, d, f, &dkb, &trlb);
    alpha = (double) (*tau - dka) / (double) (dkb - dka);
    *trl = (1 - alpha) * trla + alpha * trlb;
    return 0;
} /* ehg196_ */

/* Subroutine */ int ehg197_(long *deg, long *tau, long *d, 
	double *f, long *dk, double *trl)
{
    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    static double g1;

    *dk = 0;
    if (*deg == 1) {
	*dk = *d + 1;
    }
    if (*deg == 2) {
	*dk = (long) ((double) ((*d + 2) * (*d + 1)) / 2.);
    }
    g1 = (*d * -.08125 + .13) * *d + 1.05;
/* Computing MAX */
    d__1 = 0., d__2 = (g1 - *f) / *f;
    *trl = *dk * (max(d__1,d__2) + 1);
    return 0;
} /* ehg197_ */

/* Subroutine */ int ehg192_(double *y, long *d, long *n, long *
	nf, long *nv, long *nvmax, double *vval, double *lf, 
	long *lq)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long lq_dim1, lq_offset, lf_dim1, lf_dim2, lf_offset, vval_dim1, 
	    vval_offset, i__1, i__2, i__3;

    /* Local variables */
    static long i, j, i1, i2;
    static double i3;

    /* Parameter adjustments */
    lq_dim1 = *nvmax;
    lq_offset = lq_dim1 + 1;
    lq -= lq_offset;
    lf_dim1 = *d + 1;
    lf_dim2 = *nvmax;
    lf_offset = lf_dim1 * (lf_dim2 + 1);
    lf -= lf_offset;
    vval_dim1 = *d + 1;
    vval_offset = vval_dim1;
    vval -= vval_offset;
    --y;

    /* Function Body */
    ++execnt;
    i__1 = *nv;
    for (i2 = 1; i2 <= i__1; ++i2) {
	i__2 = *d;
	for (i1 = 0; i1 <= i__2; ++i1) {
	    vval[i1 + i2 * vval_dim1] = 0.;
/* L4: */
	}
/* L3: */
    }
    i__1 = *nv;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *nf;
	for (j = 1; j <= i__2; ++j) {
	    i3 = y[lq[i + j * lq_dim1]];
	    i__3 = *d;
	    for (i1 = 0; i1 <= i__3; ++i1) {
		vval[i1 + i * vval_dim1] += i3 * lf[i1 + (i + j * lf_dim2) * 
			lf_dim1];
/* L7: */
	    }
/* L6: */
	}
/* L5: */
    }
    return 0;
} /* ehg192_ */

double ehg128_(double *z, long *d, long *ncmax, long *vc, 
	long *a, double *xi, long *lo, long *hi, long *c, 
	double *v, long *nvmax, double *vval)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long c_dim1, c_offset, v_dim1, v_offset, vval_dim1, vval_offset, i__1, 
	    i__2;
    double ret_val, d__1;

    /* Local variables */
    extern /* Subroutine */ int ehg182_(long *), ehg184_(char *, 
	    double *, long *, long *, long);
    static double g[2304]	/* was [9][256] */, h;
    static long i, j, m;
    static double s;
    static long t[20];
    static double xibar, g0[9], g1[9];
    static long i1;
    static long i2, i3, i4, i5, i6, i7, i8, i9;
    static double v0, v1;
    static long i10;
    static long i11, i12;
    static double ge;
    static long ig, ii, lg;
    static double gn;
    static long ll;
    static double gs, gw;
    static long nt, ur;
    static double gpe, gpn, gps, gpw, sew, sns, phi0, phi1, psi0, psi1;

    /* Parameter adjustments */
    vval_dim1 = *d + 1;
    vval_offset = vval_dim1;
    vval -= vval_offset;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --hi;
    --lo;
    --xi;
    --a;
    --z;

    /* Function Body */
    ++execnt;
/*     locate enclosing cell */
    nt = 1;
    t[nt - 1] = 1;
    j = 1;
/*     top of while loop */
L3:
    if (! (a[j] != 0)) {
	goto L4;
    }
    ++nt;
    if (z[a[j]] < xi[j]) {
	i1 = lo[j];
    } else {
	i1 = hi[j];
    }
    t[nt - 1] = i1;
    if (! (nt < 20)) {
	ehg182_(&c__181);
    }
    j = t[nt - 1];
    goto L3;
/*     bottom of while loop */
L4:
/*     tensor */
    i__1 = *vc;
    for (i12 = 1; i12 <= i__1; ++i12) {
	i__2 = *d;
	for (i11 = 0; i11 <= i__2; ++i11) {
	    g[i11 + i12 * 9 - 9] = vval[i11 + c[i12 + j * c_dim1] * vval_dim1]
		    ;
/* L6: */
	}
/* L5: */
    }
    lg = *vc;
    ll = c[j * c_dim1 + 1];
    ur = c[*vc + j * c_dim1];
    for (i = *d; i >= 1; --i) {
	h = (z[i] - v[ll + i * v_dim1]) / (v[ur + i * v_dim1] - v[ll + i * 
		v_dim1]);
	if (h < -.001) {
	    ehg184_("eval ", &z[1], d, &c__1, 5L);
	    ehg184_("lowerlimit ", &v[ll + v_dim1], d, nvmax, 11L);
	} else {
	    if (1.001 < h) {
		ehg184_("eval ", &z[1], d, &c__1, 5L);
		ehg184_("upperlimit ", &v[ur + v_dim1], d, nvmax, 11L);
	    }
	}
	if (-.001 <= h) {
	    i2 = h <= 1.001;
	} else {
	    i2 = false;
	}
	if (! i2) {
	    ehg182_(&c__122);
	}
	lg = (long) ((double) lg / 2.);
	i__1 = lg;
	for (ig = 1; ig <= i__1; ++ig) {
/*           Hermite basis */
/* Computing 2nd power */
	    d__1 = 1 - h;
	    phi0 = d__1 * d__1 * (h * 2 + 1);
/* Computing 2nd power */
	    d__1 = h;
	    phi1 = d__1 * d__1 * (3 - h * 2);
/* Computing 2nd power */
	    d__1 = 1 - h;
	    psi0 = h * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = h;
	    psi1 = d__1 * d__1 * (h - 1);
	    g[ig * 9 - 9] = phi0 * g[ig * 9 - 9] + phi1 * g[(ig + lg) * 9 - 9]
		     + (psi0 * g[i + ig * 9 - 9] + psi1 * g[i + (ig + lg) * 9 
		    - 9]) * (v[ur + i * v_dim1] - v[ll + i * v_dim1]);
	    i__2 = i - 1;
	    for (ii = 1; ii <= i__2; ++ii) {
		g[ii + ig * 9 - 9] = phi0 * g[ii + ig * 9 - 9] + phi1 * g[ii 
			+ (ig + lg) * 9 - 9];
/* L9: */
	    }
/* L8: */
	}
/* L7: */
    }
    s = g[0];
/*     blending */
    if (*d == 2) {
/*        ----- North ----- */
	v0 = v[ll + v_dim1];
	v1 = v[ur + v_dim1];
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g0[i11] = vval[i11 + c[j * c_dim1 + 3] * vval_dim1];
/* L10: */
	}
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g1[i11] = vval[i11 + c[j * c_dim1 + 4] * vval_dim1];
/* L11: */
	}
	xibar = v[ur + (v_dim1 << 1)];
	m = nt - 1;
/*        top of while loop */
L12:
	if (m == 0) {
	    i4 = true;
	} else {
	    if (a[t[m - 1]] == 2) {
		i3 = xi[t[m - 1]] == xibar;
	    } else {
		i3 = false;
	    }
	    i4 = i3;
	}
	if (i4) {
	    goto L13;
	}
	--m;
/*           voidp junk */
	goto L12;
/*        bottom of while loop */
L13:
	if (m >= 1) {
	    m = hi[t[m - 1]];
/*           top of while loop */
L14:
	    if (! (a[m] != 0)) {
		goto L15;
	    }
	    if (z[a[m]] < xi[m]) {
		m = lo[m];
	    } else {
		m = hi[m];
	    }
	    goto L14;
/*           bottom of while loop */
L15:
	    if (v0 < v[c[m * c_dim1 + 1] + v_dim1]) {
		v0 = v[c[m * c_dim1 + 1] + v_dim1];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g0[i11] = vval[i11 + c[m * c_dim1 + 1] * vval_dim1];
/* L16: */
		}
	    }
	    if (v[c[m * c_dim1 + 2] + v_dim1] < v1) {
		v1 = v[c[m * c_dim1 + 2] + v_dim1];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g1[i11] = vval[i11 + c[m * c_dim1 + 2] * vval_dim1];
/* L17: */
		}
	    }
	}
	h = (z[1] - v0) / (v1 - v0);
/*        Hermite basis */
/* Computing 2nd power */
	d__1 = 1 - h;
	phi0 = d__1 * d__1 * (h * 2 + 1);
/* Computing 2nd power */
	d__1 = h;
	phi1 = d__1 * d__1 * (3 - h * 2);
/* Computing 2nd power */
	d__1 = 1 - h;
	psi0 = h * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = h;
	psi1 = d__1 * d__1 * (h - 1);
	gn = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[1] + psi1 * g1[1]) * (
		v1 - v0);
	gpn = phi0 * g0[2] + phi1 * g1[2];
/*        ----- South ----- */
	v0 = v[ll + v_dim1];
	v1 = v[ur + v_dim1];
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g0[i11] = vval[i11 + c[j * c_dim1 + 1] * vval_dim1];
/* L18: */
	}
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g1[i11] = vval[i11 + c[j * c_dim1 + 2] * vval_dim1];
/* L19: */
	}
	xibar = v[ll + (v_dim1 << 1)];
	m = nt - 1;
/*        top of while loop */
L20:
	if (m == 0) {
	    i6 = true;
	} else {
	    if (a[t[m - 1]] == 2) {
		i5 = xi[t[m - 1]] == xibar;
	    } else {
		i5 = false;
	    }
	    i6 = i5;
	}
	if (i6) {
	    goto L21;
	}
	--m;
/*           voidp junk */
	goto L20;
/*        bottom of while loop */
L21:
	if (m >= 1) {
	    m = lo[t[m - 1]];
/*           top of while loop */
L22:
	    if (! (a[m] != 0)) {
		goto L23;
	    }
	    if (z[a[m]] < xi[m]) {
		m = lo[m];
	    } else {
		m = hi[m];
	    }
	    goto L22;
/*           bottom of while loop */
L23:
	    if (v0 < v[c[m * c_dim1 + 3] + v_dim1]) {
		v0 = v[c[m * c_dim1 + 3] + v_dim1];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g0[i11] = vval[i11 + c[m * c_dim1 + 3] * vval_dim1];
/* L24: */
		}
	    }
	    if (v[c[m * c_dim1 + 4] + v_dim1] < v1) {
		v1 = v[c[m * c_dim1 + 4] + v_dim1];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g1[i11] = vval[i11 + c[m * c_dim1 + 4] * vval_dim1];
/* L25: */
		}
	    }
	}
	h = (z[1] - v0) / (v1 - v0);
/*        Hermite basis */
/* Computing 2nd power */
	d__1 = 1 - h;
	phi0 = d__1 * d__1 * (h * 2 + 1);
/* Computing 2nd power */
	d__1 = h;
	phi1 = d__1 * d__1 * (3 - h * 2);
/* Computing 2nd power */
	d__1 = 1 - h;
	psi0 = h * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = h;
	psi1 = d__1 * d__1 * (h - 1);
	gs = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[1] + psi1 * g1[1]) * (
		v1 - v0);
	gps = phi0 * g0[2] + phi1 * g1[2];
/*        ----- East ----- */
	v0 = v[ll + (v_dim1 << 1)];
	v1 = v[ur + (v_dim1 << 1)];
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g0[i11] = vval[i11 + c[j * c_dim1 + 2] * vval_dim1];
/* L26: */
	}
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g1[i11] = vval[i11 + c[j * c_dim1 + 4] * vval_dim1];
/* L27: */
	}
	xibar = v[ur + v_dim1];
	m = nt - 1;
/*        top of while loop */
L28:
	if (m == 0) {
	    i8 = true;
	} else {
	    if (a[t[m - 1]] == 1) {
		i7 = xi[t[m - 1]] == xibar;
	    } else {
		i7 = false;
	    }
	    i8 = i7;
	}
	if (i8) {
	    goto L29;
	}
	--m;
/*           voidp junk */
	goto L28;
/*        bottom of while loop */
L29:
	if (m >= 1) {
	    m = hi[t[m - 1]];
/*           top of while loop */
L30:
	    if (! (a[m] != 0)) {
		goto L31;
	    }
	    if (z[a[m]] < xi[m]) {
		m = lo[m];
	    } else {
		m = hi[m];
	    }
	    goto L30;
/*           bottom of while loop */
L31:
	    if (v0 < v[c[m * c_dim1 + 1] + (v_dim1 << 1)]) {
		v0 = v[c[m * c_dim1 + 1] + (v_dim1 << 1)];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g0[i11] = vval[i11 + c[m * c_dim1 + 1] * vval_dim1];
/* L32: */
		}
	    }
	    if (v[c[m * c_dim1 + 3] + (v_dim1 << 1)] < v1) {
		v1 = v[c[m * c_dim1 + 3] + (v_dim1 << 1)];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g1[i11] = vval[i11 + c[m * c_dim1 + 3] * vval_dim1];
/* L33: */
		}
	    }
	}
	h = (z[2] - v0) / (v1 - v0);
/*        Hermite basis */
/* Computing 2nd power */
	d__1 = 1 - h;
	phi0 = d__1 * d__1 * (h * 2 + 1);
/* Computing 2nd power */
	d__1 = h;
	phi1 = d__1 * d__1 * (3 - h * 2);
/* Computing 2nd power */
	d__1 = 1 - h;
	psi0 = h * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = h;
	psi1 = d__1 * d__1 * (h - 1);
	ge = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[2] + psi1 * g1[2]) * (
		v1 - v0);
	gpe = phi0 * g0[1] + phi1 * g1[1];
/*        ----- West ----- */
	v0 = v[ll + (v_dim1 << 1)];
	v1 = v[ur + (v_dim1 << 1)];
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g0[i11] = vval[i11 + c[j * c_dim1 + 1] * vval_dim1];
/* L34: */
	}
	i__1 = *d;
	for (i11 = 0; i11 <= i__1; ++i11) {
	    g1[i11] = vval[i11 + c[j * c_dim1 + 3] * vval_dim1];
/* L35: */
	}
	xibar = v[ll + v_dim1];
	m = nt - 1;
/*        top of while loop */
L36:
	if (m == 0) {
	    i10 = true;
	} else {
	    if (a[t[m - 1]] == 1) {
		i9 = xi[t[m - 1]] == xibar;
	    } else {
		i9 = false;
	    }
	    i10 = i9;
	}
	if (i10) {
	    goto L37;
	}
	--m;
/*           voidp junk */
	goto L36;
/*        bottom of while loop */
L37:
	if (m >= 1) {
	    m = lo[t[m - 1]];
/*           top of while loop */
L38:
	    if (! (a[m] != 0)) {
		goto L39;
	    }
	    if (z[a[m]] < xi[m]) {
		m = lo[m];
	    } else {
		m = hi[m];
	    }
	    goto L38;
/*           bottom of while loop */
L39:
	    if (v0 < v[c[m * c_dim1 + 2] + (v_dim1 << 1)]) {
		v0 = v[c[m * c_dim1 + 2] + (v_dim1 << 1)];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g0[i11] = vval[i11 + c[m * c_dim1 + 2] * vval_dim1];
/* L40: */
		}
	    }
	    if (v[c[m * c_dim1 + 4] + (v_dim1 << 1)] < v1) {
		v1 = v[c[m * c_dim1 + 4] + (v_dim1 << 1)];
		i__1 = *d;
		for (i11 = 0; i11 <= i__1; ++i11) {
		    g1[i11] = vval[i11 + c[m * c_dim1 + 4] * vval_dim1];
/* L41: */
		}
	    }
	}
	h = (z[2] - v0) / (v1 - v0);
/*        Hermite basis */
/* Computing 2nd power */
	d__1 = 1 - h;
	phi0 = d__1 * d__1 * (h * 2 + 1);
/* Computing 2nd power */
	d__1 = h;
	phi1 = d__1 * d__1 * (3 - h * 2);
/* Computing 2nd power */
	d__1 = 1 - h;
	psi0 = h * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = h;
	psi1 = d__1 * d__1 * (h - 1);
	gw = phi0 * g0[0] + phi1 * g1[0] + (psi0 * g0[2] + psi1 * g1[2]) * (
		v1 - v0);
	gpw = phi0 * g0[1] + phi1 * g1[1];
/*        NS */
	h = (z[2] - v[ll + (v_dim1 << 1)]) / (v[ur + (v_dim1 << 1)] - v[ll + (
		v_dim1 << 1)]);
/*        Hermite basis */
/* Computing 2nd power */
	d__1 = 1 - h;
	phi0 = d__1 * d__1 * (h * 2 + 1);
/* Computing 2nd power */
	d__1 = h;
	phi1 = d__1 * d__1 * (3 - h * 2);
/* Computing 2nd power */
	d__1 = 1 - h;
	psi0 = h * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = h;
	psi1 = d__1 * d__1 * (h - 1);
	sns = phi0 * gs + phi1 * gn + (psi0 * gps + psi1 * gpn) * (v[ur + (
		v_dim1 << 1)] - v[ll + (v_dim1 << 1)]);
/*        EW */
	h = (z[1] - v[ll + v_dim1]) / (v[ur + v_dim1] - v[ll + v_dim1]);
/*        Hermite basis */
/* Computing 2nd power */
	d__1 = 1 - h;
	phi0 = d__1 * d__1 * (h * 2 + 1);
/* Computing 2nd power */
	d__1 = h;
	phi1 = d__1 * d__1 * (3 - h * 2);
/* Computing 2nd power */
	d__1 = 1 - h;
	psi0 = h * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = h;
	psi1 = d__1 * d__1 * (h - 1);
	sew = phi0 * gw + phi1 * ge + (psi0 * gpw + psi1 * gpe) * (v[ur + 
		v_dim1] - v[ll + v_dim1]);
	s = sns + sew - s;
    }
    ret_val = s;
    return ret_val;
} /* ehg128_ */

long ifloor_(double *x)
{
    /* System generated locals */
    long ret_val;

    ret_val = (long) (*x);
    if ((double) ret_val > *x) {
	--ret_val;
    }
    return ret_val;
} /* ifloor_ */

double dsign_(double *a1, double *a2)
{
    /* System generated locals */
    double ret_val;

    ret_val = abs(*a1);
    if (*a2 >= 0.) {
	ret_val = -ret_val;
    }
    return ret_val;
} /* dsign_ */

/* Subroutine */ int ehg136_(double *u, long *lm, long *m, long *
	n, long *d, long *nf, double *f, double *x, long *
	psi, double *y, double *rw, long *kernel, long *k, 
	double *dist, double *eta, double *b, long *od, 
	double *o, long *ihat, double *w, double *rcond, 
	long *sing, long *dd, long *tdeg, long *cdeg, double *
	s)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long o_dim1, o_offset, b_dim1, b_offset, s_dim1, s_offset, u_dim1, 
	    u_offset, x_dim1, x_offset, i__1, i__2, i__3;

    /* Local variables */
    extern double ddot_(long *, double *, long *, double *, 
	    long *);
    static long info;
    static double work[15];
    extern /* Subroutine */ int ehg127_(double *, long *, long *, 
	    long *, double *, double *, long *, double *, 
	    double *, long *, long *, double *, double *, 
	    double *, long *, double *, double *, long *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, long *, long *, 
	    long *, double *), ehg182_(long *);
    static double e[225]	/* was [15][15] */, g[225]	/* was [15][
	    15] */;
    static long i, j, l;
    static double q[8], scale, sigma[15];
    extern /* Subroutine */ int dqrsl_(double *, long *, long *, 
	    long *, double *, double *, double *, double *,
	     double *, double *, double *, long *, long *);
    static long i1;
    static double i2, qraux[15], dgamma[15];
    static long identi;
    static double tol;

    /* Parameter adjustments */
    s_dim1 = *od + 1;
    s_offset = s_dim1;
    s -= s_offset;
    --cdeg;
    --w;
    o_dim1 = *m;
    o_offset = o_dim1 + 1;
    o -= o_offset;
    b_dim1 = *nf;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    --eta;
    --dist;
    --rw;
    --y;
    --psi;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    u_dim1 = *lm;
    u_offset = u_dim1 + 1;
    u -= u_offset;

    /* Function Body */
/*     V -> g */
/*     U -> e */
/*     Identity -> identi */
/*     L -> o */
/*     X -> b */
    ++execnt;
    if (! (*k <= *nf - 1)) {
	ehg182_(&c__104);
    }
    if (! (*k <= 15)) {
	ehg182_(&c__105);
    }
    i__1 = *n;
    for (identi = 1; identi <= i__1; ++identi) {
	psi[identi] = identi;
/* L3: */
    }
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *d;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    q[i1 - 1] = u[l + i1 * u_dim1];
/* L5: */
	}
	ehg127_(q, n, d, nf, f, &x[x_offset], &psi[1], &y[1], &rw[1], kernel, 
		k, &dist[1], &eta[1], &b[b_offset], od, &w[1], rcond, sing, 
		sigma, e, g, dgamma, qraux, work, &tol, dd, tdeg, &cdeg[1], &
		s[l * s_dim1]);
	if (*ihat == 1) {
/*           $L sub {l,l} = */
/*           V sub {1,:} SIGMA sup {+} U sup T */
/*           (Q sup T W e sub i )$ */
	    if (! (*m == *n)) {
		ehg182_(&c__123);
	    }
/*           find $i$ such that $l = psi sub i$ */
	    i = 1;
/*           top of while loop */
L6:
	    if (! (l != psi[i])) {
		goto L7;
	    }
	    ++i;
	    if (! (i < *nf)) {
		ehg182_(&c__123);
	    }
	    goto L6;
/*           bottom of while loop */
L7:
	    i__2 = *nf;
	    for (i1 = 1; i1 <= i__2; ++i1) {
		eta[i1] = 0.;
/* L8: */
	    }
	    eta[i] = w[i];
/*           $eta = Q sup T W e sub i$ */
	    dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], &eta[1], &eta[1], 
		    &eta[1], &eta[1], &eta[1], &c__1000, &info);
/*           $gamma = U sup T eta sub {1:k}$ */
	    i__2 = *k;
	    for (i1 = 1; i1 <= i__2; ++i1) {
		dgamma[i1 - 1] = 0.;
/* L9: */
	    }
	    i__2 = *k;
	    for (j = 1; j <= i__2; ++j) {
		i2 = eta[j];
		i__3 = *k;
		for (i1 = 1; i1 <= i__3; ++i1) {
		    dgamma[i1 - 1] += i2 * e[j + i1 * 15 - 16];
/* L11: */
		}
/* L10: */
	    }
/*           $gamma = SIGMA sup {+} gamma$ */
	    i__2 = *k;
	    for (j = 1; j <= i__2; ++j) {
		if (tol < sigma[j - 1]) {
		    dgamma[j - 1] /= sigma[j - 1];
		} else {
		    dgamma[j - 1] = 0.;
		}
/* L12: */
	    }
/*           voidp junk */
/*           voidp junk */
	    o[l + o_dim1] = ddot_(k, g, &c__15, dgamma, &c__1);
	} else {
	    if (*ihat == 2) {
/*              $L sub {l,:} = */
/*              V sub {1,:} SIGMA sup {+} */
/*              ( U sup T Q sup T ) W $ */
		i__2 = *n;
		for (i1 = 1; i1 <= i__2; ++i1) {
		    o[l + i1 * o_dim1] = 0.;
/* L13: */
		}
		i__2 = *k;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = *nf;
		    for (i1 = 1; i1 <= i__3; ++i1) {
			eta[i1] = 0.;
/* L15: */
		    }
		    i__3 = *k;
		    for (i1 = 1; i1 <= i__3; ++i1) {
			eta[i1] = e[i1 + j * 15 - 16];
/* L16: */
		    }
		    dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], &eta[1], 
			    work, work, work, work, &c__10000, &info);
		    if (tol < sigma[j - 1]) {
			scale = 1. / sigma[j - 1];
		    } else {
			scale = 0.;
		    }
		    i__3 = *nf;
		    for (i1 = 1; i1 <= i__3; ++i1) {
			eta[i1] *= scale * w[i1];
/* L17: */
		    }
		    i__3 = *nf;
		    for (i = 1; i <= i__3; ++i) {
			o[l + psi[i] * o_dim1] += g[j * 15 - 15] * eta[i];
/* L18: */
		    }
/* L14: */
		}
	    }
	}
/* L4: */
    }
    return 0;
} /* ehg136_ */


/* Subroutine */ int ehg139_(double *v, long *nvmax, long *nv, 
	long *n, long *d, long *nf, double *f, double *x, 
	long *pi, long *psi, double *y, double *rw, double *
	trl, long *kernel, long *k, double *dist, long *phi, 
	double *eta, double *b, long *od, double *w, 
	double *diagl, double *vval2, long *ncmax, long *vc, 
	long *a, double *xi, long *lo, long *hi, long *c, 
	long *vhit, double *rcond, long *sing, long *dd, long 
	*tdeg, long *cdeg, long *lq, double *lf, long *setlf, 
	double *s)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long lq_dim1, lq_offset, c_dim1, c_offset, lf_dim1, lf_dim2, lf_offset,
	     b_dim1, b_offset, s_dim1, s_offset, v_dim1, v_offset, vval2_dim1,
	     vval2_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static long leaf[256];
    extern double ddot_(long *, double *, long *, double *, 
	    long *);
    static long info;
    static double term, work[15];
    extern /* Subroutine */ int ehg127_(double *, long *, long *, 
	    long *, double *, double *, long *, double *, 
	    double *, long *, long *, double *, double *, 
	    double *, long *, double *, double *, long *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, long *, long *, 
	    long *, double *), ehg182_(long *), ehg137_(double *
	    , long *, long *, long *, long *, long *, long *
	    , long *, long *, long *, double *, long *, 
	    long *, long *, double *);
    extern double ehg128_(double *, long *, long *, long *, 
	    long *, double *, long *, long *, long *, 
	    double *, long *, double *);
    static double e[225]	/* was [15][15] */;
    static long i, j, l, ileaf;
    static double q[8];
    static long nleaf;
    static double scale, u[225]	/* was [15][15] */, z[8], sigma[15];
    extern /* Subroutine */ int dqrsl_(double *, long *, long *, 
	    long *, double *, double *, double *, double *,
	     double *, double *, double *, long *, long *);
    static double i1;
    static long i2, i3;
    static double i4;
    static long i5, i6;
    static double i7, qraux[15];
    static long ii;
    static double dgamma[15];
    static long identi;
    static double tol;

    /* Parameter adjustments */
    s_dim1 = *od + 1;
    s_offset = s_dim1;
    s -= s_offset;
    lf_dim1 = *d + 1;
    lf_dim2 = *nvmax;
    lf_offset = lf_dim1 * (lf_dim2 + 1);
    lf -= lf_offset;
    lq_dim1 = *nvmax;
    lq_offset = lq_dim1 + 1;
    lq -= lq_offset;
    --cdeg;
    --vhit;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --hi;
    --lo;
    --xi;
    --a;
    vval2_dim1 = *d + 1;
    vval2_offset = vval2_dim1;
    vval2 -= vval2_offset;
    --diagl;
    --w;
    b_dim1 = *nf;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    --eta;
    --phi;
    --dist;
    --rw;
    --y;
    --psi;
    --pi;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;

    /* Function Body */
/*     V -> e */
/*     Identity -> identi */
/*     X -> b */
    ++execnt;
/*     l2fit with trace(L) */
    if (! (*k <= *nf - 1)) {
	ehg182_(&c__104);
    }
    if (! (*k <= 15)) {
	ehg182_(&c__105);
    }
    if (*trl != 0.) {
	i__1 = *n;
	for (i5 = 1; i5 <= i__1; ++i5) {
	    diagl[i5] = 0.;
/* L3: */
	}
	i__1 = *nv;
	for (i6 = 1; i6 <= i__1; ++i6) {
	    i__2 = *d;
	    for (i5 = 0; i5 <= i__2; ++i5) {
		vval2[i5 + i6 * vval2_dim1] = 0.;
/* L5: */
	    }
/* L4: */
	}
    }
    i__1 = *n;
    for (identi = 1; identi <= i__1; ++identi) {
	psi[identi] = identi;
/* L6: */
    }
    i__1 = *nv;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *d;
	for (i5 = 1; i5 <= i__2; ++i5) {
	    q[i5 - 1] = v[l + i5 * v_dim1];
/* L8: */
	}
	ehg127_(q, n, d, nf, f, &x[x_offset], &psi[1], &y[1], &rw[1], kernel, 
		k, &dist[1], &eta[1], &b[b_offset], od, &w[1], rcond, sing, 
		sigma, u, e, dgamma, qraux, work, &tol, dd, tdeg, &cdeg[1], &
		s[l * s_dim1]);
	if (*trl != 0.) {
/*           invert $psi$ */
	    i__2 = *n;
	    for (i5 = 1; i5 <= i__2; ++i5) {
		phi[i5] = 0;
/* L9: */
	    }
	    i__2 = *nf;
	    for (i = 1; i <= i__2; ++i) {
		phi[psi[i]] = i;
/* L10: */
	    }
	    i__2 = *d;
	    for (i5 = 1; i5 <= i__2; ++i5) {
		z[i5 - 1] = v[l + i5 * v_dim1];
/* L11: */
	    }
	    ehg137_(z, &vhit[l], leaf, &nleaf, d, nv, nvmax, ncmax, vc, &a[1],
		     &xi[1], &lo[1], &hi[1], &c[c_offset], &v[v_offset]);
	    i__2 = nleaf;
	    for (ileaf = 1; ileaf <= i__2; ++ileaf) {
		i__3 = hi[leaf[ileaf - 1]];
		for (ii = lo[leaf[ileaf - 1]]; ii <= i__3; ++ii) {
		    i = phi[pi[ii]];
		    if (i != 0) {
			if (! (psi[i] == pi[ii])) {
			    ehg182_(&c__194);
			}
			i__4 = *nf;
			for (i5 = 1; i5 <= i__4; ++i5) {
			    eta[i5] = 0.;
/* L14: */
			}
			eta[i] = w[i];
/*                    $eta = Q sup T W e sub i$ */
			dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], work, 
				&eta[1], &eta[1], work, work, &c__1000, &info)
				;
			i__4 = *k;
			for (j = 1; j <= i__4; ++j) {
			    if (tol < sigma[j - 1]) {
				i4 = ddot_(k, &u[j * 15 - 15], &c__1, &eta[1],
					 &c__1) / sigma[j - 1];
			    } else {
				i4 = 0.;
			    }
			    dgamma[j - 1] = i4;
/* L15: */
			}
			i__4 = *d + 1;
			for (j = 1; j <= i__4; ++j) {
			    vval2[j - 1 + l * vval2_dim1] = ddot_(k, &e[j - 1]
				    , &c__15, dgamma, &c__1);
/* L16: */
			}
			i__4 = *d;
			for (i5 = 1; i5 <= i__4; ++i5) {
			    z[i5 - 1] = x[pi[ii] + i5 * x_dim1];
/* L17: */
			}
			term = ehg128_(z, d, ncmax, vc, &a[1], &xi[1], &lo[1],
				 &hi[1], &c[c_offset], &v[v_offset], nvmax, &
				vval2[vval2_offset]);
			diagl[pi[ii]] += term;
			i__4 = *d;
			for (i5 = 0; i5 <= i__4; ++i5) {
			    vval2[i5 + l * vval2_dim1] = 0.;
/* L18: */
			}
		    }
/* L13: */
		}
/* L12: */
	    }
	}
	if (*setlf) {
/*           $Lf sub {:,l,:} = V SIGMA sup {+} U sup T Q sup T W$ 
*/
	    if (! (*k >= *d + 1)) {
		ehg182_(&c__196);
	    }
	    i__2 = *nf;
	    for (i5 = 1; i5 <= i__2; ++i5) {
		lq[l + i5 * lq_dim1] = psi[i5];
/* L19: */
	    }
	    i__2 = *nf;
	    for (i6 = 1; i6 <= i__2; ++i6) {
		i__3 = *d;
		for (i5 = 0; i5 <= i__3; ++i5) {
		    lf[i5 + (l + i6 * lf_dim2) * lf_dim1] = 0.;
/* L21: */
		}
/* L20: */
	    }
	    i__2 = *k;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *nf;
		for (i5 = 1; i5 <= i__3; ++i5) {
		    eta[i5] = 0.;
/* L23: */
		}
		i__3 = *k;
		for (i5 = 1; i5 <= i__3; ++i5) {
		    eta[i5] = u[i5 + j * 15 - 16];
/* L24: */
		}
		dqrsl_(&b[b_offset], nf, nf, k, qraux, &eta[1], &eta[1], work,
			 work, work, work, &c__10000, &info);
		if (tol < sigma[j - 1]) {
		    scale = 1. / sigma[j - 1];
		} else {
		    scale = 0.;
		}
		i__3 = *nf;
		for (i5 = 1; i5 <= i__3; ++i5) {
		    eta[i5] *= scale * w[i5];
/* L25: */
		}
		i__3 = *nf;
		for (i = 1; i <= i__3; ++i) {
		    i7 = eta[i];
		    i__4 = *d;
		    for (i5 = 0; i5 <= i__4; ++i5) {
			lf[i5 + (l + i * lf_dim2) * lf_dim1] += e[i5 + 1 + j *
				 15 - 16] * i7;
/* L27: */
		    }
/* L26: */
		}
/* L22: */
	    }
	}
/* L7: */
    }
    if (*trl != 0.) {
	if (*n <= 0) {
	    *trl = 0.;
	} else {
	    i3 = *n;
	    i1 = diagl[i3];
	    for (i2 = i3 - 1; i2 >= 1; --i2) {
		i1 = diagl[i2] + i1;
/* L28: */
	    }
	    *trl = i1;
	}
    }
    return 0;
} /* ehg139_ */

/* Subroutine */ int dqrdc_(double *x, long *ldx, long *n, long *
	p, double *qraux, long *jpvt, double *work, long *job)
{
    /* System generated locals */
    long x_dim1, x_offset, i__1, i__2, i__3;
    double d__1, d__2;

    /* Builtin functions */
    double d_sign(double *, double *), sqrt(double);

    /* Local variables */
    static long negj;
    extern double ddot_(long *, double *, long *, double *, 
	    long *);
    static long maxj;
    extern double dnrm2_(long *, double *, long *);
    static long j, l;
    static double t;
    extern /* Subroutine */ int dscal_(long *, double *, double *, 
	    long *), dswap_(long *, double *, long *, double 
	    *, long *);
    static long swapj;
    extern /* Subroutine */ int daxpy_(long *, double *, double *, 
	    long *, double *, long *);
    static double nrmxl;
    static long jj, jp, pl, pu;
    static double tt, maxnrm;
    static long lp1, lup;


/*     dqrdc uses householder transformations to compute the qr */
/*     factorization of an n by p matrix x.  column pivoting */
/*     based on the 2-norms of the reduced columns may be */
/*     performed at the users option. */

/*     on entry */

/*        x       double precision(ldx,p), where ldx .ge. n. */
/*                x contains the matrix whose decomposition is to be */
/*                computed. */

/*        ldx     long. */
/*                ldx is the leading dimension of the array x. */

/*        n       long. */
/*                n is the number of rows of the matrix x. */

/*        p       long. */
/*                p is the number of columns of the matrix x. */

/*        jpvt    long(p). */
/*                jpvt contains longs that control the selection */
/*                of the pivot columns.  the k-th column x(k) of x */
/*                is placed in one of three classes according to the */
/*                value of jpvt(k). */

/*                   if jpvt(k) .gt. 0, then x(k) is an initial */
/*                                      column. */

/*                   if jpvt(k) .eq. 0, then x(k) is a free column. */

/*                   if jpvt(k) .lt. 0, then x(k) is a final column. */

/*                before the decomposition is computed, initial columns */
/*                are moved to the beginning of the array x and final */
/*                columns to the end.  both initial and final columns */
/*                are frozen in place during the computation and only */
/*                free columns are moved.  at the k-th stage of the */
/*                reduction, if x(k) is occupied by a free column */
/*                it is interchanged with the free column of largest */
/*                reduced norm.  jpvt is not referenced if */
/*                job .eq. 0. */

/*        work    double precision(p). */
/*                work is a work array.  work is not referenced if */
/*                job .eq. 0. */

/*        job     long. */
/*                job is an long that initiates column pivoting. */
/*                if job .eq. 0, no pivoting is done. */
/*                if job .ne. 0, pivoting is done. */

/*     on return */

/*        x       x contains in its upper triangle the upper */
/*                triangular matrix r of the qr factorization. */
/*                below its diagonal x contains information from */
/*                which the orthogonal part of the decomposition */
/*                can be recovered.  note that if pivoting has */
/*                been requested, the decomposition is not that */
/*                of the original matrix x but that of x */
/*                with its columns permuted as described by jpvt. */

/*        qraux   double precision(p). */
/*                qraux contains further information required to recover 
*/
/*                the orthogonal part of the decomposition. */

/*        jpvt    jpvt(k) contains the index of the column of the */
/*                original matrix that has been interchanged into */
/*                the k-th column, if pivoting was requested. */

/*     linpack. this version dated 08/14/78 . */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     dqrdc uses the following functions and subprograms. */

/*     blas daxpy,ddot,dscal,dswap,dnrm2 */
/*     fortran dabs,dmax1,min0,dsqrt */

/*     internal variables */



    /* Parameter adjustments */
    --work;
    --jpvt;
    --qraux;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    pl = 1;
    pu = 0;
    if (*job == 0) {
	goto L60;
    }

/*        pivoting has been requested.  rearrange the columns */
/*        according to jpvt. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	swapj = jpvt[j] > 0;
	negj = jpvt[j] < 0;
	jpvt[j] = j;
	if (negj) {
	    jpvt[j] = -j;
	}
	if (! swapj) {
	    goto L10;
	}
	if (j != pl) {
	    dswap_(n, &x[pl * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	}
	jpvt[j] = jpvt[pl];
	jpvt[pl] = j;
	++pl;
L10:
/* L20: */
	;
    }
    pu = *p;
    i__1 = *p;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *p - jj + 1;
	if (jpvt[j] >= 0) {
	    goto L40;
	}
	jpvt[j] = -jpvt[j];
	if (j == pu) {
	    goto L30;
	}
	dswap_(n, &x[pu * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	jp = jpvt[pu];
	jpvt[pu] = jpvt[j];
	jpvt[j] = jp;
L30:
	--pu;
L40:
/* L50: */
	;
    }
L60:

/*     compute the norms of the free columns. */

    if (pu < pl) {
	goto L80;
    }
    i__1 = pu;
    for (j = pl; j <= i__1; ++j) {
	qraux[j] = dnrm2_(n, &x[j * x_dim1 + 1], &c__1);
	work[j] = qraux[j];
/* L70: */
    }
L80:

/*     perform the householder reduction of x. */

    lup = min(*n,*p);
    i__1 = lup;
    for (l = 1; l <= i__1; ++l) {
	if (l < pl || l >= pu) {
	    goto L120;
	}

/*           locate the column of largest norm and bring it */
/*           into the pivot position. */

	maxnrm = 0.;
	maxj = l;
	i__2 = pu;
	for (j = l; j <= i__2; ++j) {
	    if (qraux[j] <= maxnrm) {
		goto L90;
	    }
	    maxnrm = qraux[j];
	    maxj = j;
L90:
/* L100: */
	    ;
	}
	if (maxj == l) {
	    goto L110;
	}
	dswap_(n, &x[l * x_dim1 + 1], &c__1, &x[maxj * x_dim1 + 1], &c__1);
	qraux[maxj] = qraux[l];
	work[maxj] = work[l];
	jp = jpvt[maxj];
	jpvt[maxj] = jpvt[l];
	jpvt[l] = jp;
L110:
L120:
	qraux[l] = 0.;
	if (l == *n) {
	    goto L190;
	}

/*           compute the householder transformation for column l. */

	i__2 = *n - l + 1;
	nrmxl = dnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (nrmxl == 0.) {
	    goto L180;
	}
	if (x[l + l * x_dim1] != 0.) {
	    nrmxl = d_sign(&nrmxl, &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / nrmxl;
	dscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;

/*              apply the transformation to the remaining columns, */
/*              updating the norms. */

	lp1 = l + 1;
	if (*p < lp1) {
	    goto L170;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
	    if (j < pl || j > pu) {
		goto L150;
	    }
	    if (qraux[j] == 0.) {
		goto L150;
	    }
/* Computing 2nd power */
	    d__2 = (d__1 = x[l + j * x_dim1], abs(d__1)) / qraux[j];
	    tt = 1. - d__2 * d__2;
	    tt = max(tt,0.);
	    t = tt;
/* Computing 2nd power */
	    d__1 = qraux[j] / work[j];
	    tt = tt * .05 * (d__1 * d__1) + 1.;
	    if (tt == 1.) {
		goto L130;
	    }
	    qraux[j] *= sqrt(t);
	    goto L140;
L130:
	    i__3 = *n - l;
	    qraux[j] = dnrm2_(&i__3, &x[l + 1 + j * x_dim1], &c__1);
	    work[j] = qraux[j];
L140:
L150:
/* L160: */
	    ;
	}
L170:

/*              save the transformation. */

	qraux[l] = x[l + l * x_dim1];
	x[l + l * x_dim1] = -nrmxl;
L180:
L190:
/* L200: */
	;
    }
    return 0;
} /* dqrdc_ */

long idamax_(long *n, double *dx, long *incx)
{
    /* System generated locals */
    long ret_val, i__1;
    double d__1;

    /* Local variables */
    static double dmax_;
    static long i, ix;


/*     finds the index of element having max. absolute value. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    dmax_ = abs(dx[1]);
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((d__1 = dx[ix], abs(d__1)) <= dmax_) {
	    goto L5;
	}
	ret_val = i;
	dmax_ = (d__1 = dx[ix], abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    dmax_ = abs(dx[1]);
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if ((d__1 = dx[i], abs(d__1)) <= dmax_) {
	    goto L30;
	}
	ret_val = i;
	dmax_ = (d__1 = dx[i], abs(d__1));
L30:
	;
    }
    return ret_val;
} /* idamax_ */

/* Subroutine */ int lowesb(double *xx, double *yy, double *ww, 
	double *diagl, long *infl, long *iv, long *liv, long *
	lv, double *wv)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long i__1;
    double d__1;

    /* Local variables */
    extern /* Subroutine */ int ehg131_(double *, double *, 
	    double *, double *, double *, long *, long *, 
	    long *, long *, long *, long *, long *, long *, 
	    long *, long *, double *, long *, long *, long 
	    *, long *, long *, long *, double *, long *, 
	    double *, double *, double *, double *, 
	    double *, long *, double *, double *, double *,
	     double *, long *, long *, long *, long *, 
	    long *, double *, long *), ehg182_(long *), ehg183_(
	    char *, long *, long *, long *, long);
    static long setlf;
    extern long ifloor_(double *);
    static double trl;

    /* Parameter adjustments */
    --wv;
    --iv;
    --diagl;
    --ww;
    --yy;
    --xx;

    /* Function Body */
    ++execnt;
    if (! (iv[28] != 173)) {
	ehg182_(&c__174);
    }
    if (iv[28] != 172) {
	if (! (iv[28] == 171)) {
	    ehg182_(&c__171);
	}
    }
    iv[28] = 173;
    if (*infl) {
	trl = 1.;
    } else {
	trl = 0.;
    }
    setlf = iv[27] != iv[25];
    d__1 = iv[3] * wv[2];
    i__1 = ifloor_(&d__1);
    ehg131_(&xx[1], &yy[1], &ww[1], &trl, &diagl[1], &iv[20], &iv[29], &iv[3],
	     &iv[2], &iv[5], &iv[17], &iv[4], &iv[6], &iv[14], &iv[19], &wv[1]
	    , &iv[iv[7]], &iv[iv[8]], &iv[iv[9]], &iv[iv[10]], &iv[iv[22]], &
	    iv[iv[27]], &wv[iv[11]], &iv[iv[23]], &wv[iv[13]], &wv[iv[12]], &
	    wv[iv[15]], &wv[iv[16]], &wv[iv[18]], &i__1, &wv[3], &wv[iv[26]], 
	    &wv[iv[24]], &wv[4], &iv[30], &iv[33], &iv[32], &iv[41], &iv[iv[
	    25]], &wv[iv[34]], &setlf);
    if ((double) iv[14] < iv[6] + (double) iv[4] / 2.) {
	ehg183_("Warning. k-d tree limited by memory; nvmax=", &iv[14], &c__1,
		 &c__1, 43L);
    } else {
	if (iv[17] < iv[5] + 2) {
	    ehg183_("Warning. k-d tree limited by memory. ncmax=", &iv[17], &
		    c__1, &c__1, 43L);
	}
    }
    return 0;
} /* lowesb_ */

/* Subroutine */ int lowesd(long *versio, long *iv, long *liv, 
	long *lv, double *v, long *d, long *n, double *f, 
	long *ideg, long *nvmax, long *setlf)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long i__1, i__2;
    double d__1;

    /* Builtin functions */
    long pow_ii(long *, long *);

    /* Local variables */
    extern /* Subroutine */ int ehg182_(long *);
    static long i, j, ncmax, bound, i1, i2, nf, vc;
    extern long ifloor_(double *);

    /* Parameter adjustments */
    --v;
    --iv;

    /* Function Body */
/*     version -> versio */
    ++execnt;
    if (! (*versio == 106)) {
	ehg182_(&c__100);
    }
    iv[28] = 171;
    iv[2] = *d;
    iv[3] = *n;
    vc = pow_ii(&c__2, d);
    iv[4] = vc;
    if (! (0. < *f)) {
	ehg182_(&c__120);
    }
/* Computing MIN */
    d__1 = *n * *f;
    i__1 = *n, i__2 = ifloor_(&d__1);
    nf = min(i__1,i__2);
    iv[19] = nf;
    iv[20] = 1;
    if (*ideg == 0) {
	i1 = 1;
    } else {
	if (*ideg == 1) {
	    i1 = *d + 1;
	} else {
	    if (*ideg == 2) {
		i1 = (long) ((double) ((*d + 2) * (*d + 1)) / 2.);
	    }
	}
    }
    iv[29] = i1;
    iv[21] = 1;
    iv[14] = *nvmax;
    ncmax = *nvmax;
    iv[17] = ncmax;
    iv[30] = 0;
    iv[32] = *ideg;
    if (! (*ideg >= 0)) {
	ehg182_(&c__195);
    }
    if (! (*ideg <= 2)) {
	ehg182_(&c__195);
    }
    iv[33] = *d;
    for (i2 = 41; i2 <= 49; ++i2) {
	iv[i2] = *ideg;
/* L3: */
    }
    iv[7] = 50;
    iv[8] = iv[7] + ncmax;
    iv[9] = iv[8] + vc * ncmax;
    iv[10] = iv[9] + ncmax;
    iv[22] = iv[10] + ncmax;
/*     initialize permutation */
    j = iv[22] - 1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	iv[j + i] = i;
/* L4: */
    }
    iv[23] = iv[22] + *n;
    iv[25] = iv[23] + *nvmax;
    if (*setlf) {
	iv[27] = iv[25] + *nvmax * nf;
    } else {
	iv[27] = iv[25];
    }
    bound = iv[27] + *n;
    if (! (bound - 1 <= *liv)) {
	ehg182_(&c__102);
    }
    iv[11] = 50;
    iv[13] = iv[11] + *nvmax * *d;
    iv[12] = iv[13] + (*d + 1) * *nvmax;
    iv[15] = iv[12] + ncmax;
    iv[16] = iv[15] + *n;
    iv[18] = iv[16] + nf;
    iv[24] = iv[18] + iv[29] * nf;
    iv[34] = iv[24] + (*d + 1) * *nvmax;
    if (*setlf) {
	iv[26] = iv[34] + (*d + 1) * *nvmax * nf;
    } else {
	iv[26] = iv[34];
    }
    bound = iv[26] + nf;
    if (! (bound - 1 <= *lv)) {
	ehg182_(&c__103);
    }
    v[1] = *f;
    v[2] = .05;
    v[3] = 0.;
    v[4] = 1.;
    return 0;
} /* lowesd_ */

/* Subroutine */ int lowese(long *iv, long *liv, long *lv, 
	double *wv, long *m, double *z, double *s)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long z_dim1, z_offset;

    /* Local variables */
    extern /* Subroutine */ int ehg133_(long *, long *, long *, 
	    long *, long *, long *, long *, long *, long *, 
	    long *, double *, double *, double *, long *, 
	    double *, double *), ehg182_(long *);

    /* Parameter adjustments */
    --s;
    z_dim1 = *m;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --wv;
    --iv;

    /* Function Body */
    ++execnt;
    if (! (iv[28] != 172)) {
	ehg182_(&c__172);
    }
    if (! (iv[28] == 173)) {
	ehg182_(&c__173);
    }
    ehg133_(&iv[3], &iv[2], &iv[4], &iv[14], &iv[5], &iv[17], &iv[iv[7]], &iv[
	    iv[8]], &iv[iv[9]], &iv[iv[10]], &wv[iv[11]], &wv[iv[13]], &wv[iv[
	    12]], m, &z[z_offset], &s[1]);
    return 0;
} /* lowese_ */

/* Subroutine */ int lowesf(double *xx, double *yy, double *ww, 
	long *iv, long *liv, long *lv, double *wv, long *m, 
	double *z, double *l, long *ihat, double *s)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long l_dim1, l_offset, z_dim1, z_offset;

    /* Local variables */
    extern /* Subroutine */ int ehg136_(double *, long *, long *, 
	    long *, long *, long *, double *, double *, 
	    long *, double *, double *, long *, long *, 
	    double *, double *, double *, long *, double *,
	     long *, double *, double *, long *, long *, 
	    long *, long *, double *), ehg182_(long *);
    static long i1;

    /* Parameter adjustments */
    --s;
    l_dim1 = *m;
    l_offset = l_dim1 + 1;
    l -= l_offset;
    z_dim1 = *m;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --wv;
    --iv;
    --ww;
    --yy;
    --xx;

    /* Function Body */
    ++execnt;
    if (171 <= iv[28]) {
	i1 = iv[28] <= 174;
    } else {
	i1 = false;
    }
    if (! i1) {
	ehg182_(&c__171);
    }
    iv[28] = 172;
    if (! (iv[14] >= iv[19])) {
	ehg182_(&c__186);
    }
    ehg136_(&z[z_offset], m, m, &iv[3], &iv[2], &iv[19], &wv[1], &xx[1], &iv[
	    iv[22]], &yy[1], &ww[1], &iv[20], &iv[29], &wv[iv[15]], &wv[iv[16]
	    ], &wv[iv[18]], &c__0, &l[l_offset], ihat, &wv[iv[26]], &wv[4], &
	    iv[30], &iv[33], &iv[32], &iv[41], &s[1]);
    return 0;
} /* lowesf_ */

/* Subroutine */ int lowesl(long *iv, long *liv, long *lv, 
	double *wv, long *m, double *z, double *l)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long l_dim1, l_offset, z_dim1, z_offset;

    /* Local variables */
    extern /* Subroutine */ int ehg182_(long *), ehg191_(long *, 
	    double *, double *, long *, long *, long *, 
	    long *, long *, long *, long *, double *, long 
	    *, long *, long *, double *, long *, double *, 
	    double *, long *);

    /* Parameter adjustments */
    l_dim1 = *m;
    l_offset = l_dim1 + 1;
    l -= l_offset;
    z_dim1 = *m;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --wv;
    --iv;

    /* Function Body */
    ++execnt;
    if (! (iv[28] != 172)) {
	ehg182_(&c__172);
    }
    if (! (iv[28] == 173)) {
	ehg182_(&c__173);
    }
    if (! (iv[26] != iv[34])) {
	ehg182_(&c__175);
    }
    ehg191_(m, &z[z_offset], &l[l_offset], &iv[2], &iv[3], &iv[19], &iv[6], &
	    iv[17], &iv[4], &iv[iv[7]], &wv[iv[12]], &iv[iv[10]], &iv[iv[9]], 
	    &iv[iv[8]], &wv[iv[11]], &iv[14], &wv[iv[24]], &wv[iv[34]], &iv[
	    iv[25]]);
    return 0;
} /* lowesl_ */

/* Subroutine */ int lowesr_(double *yy, long *iv, long *liv, 
	long *lv, double *wv)
{
    /* Initialized data */

    static long execnt = 0;

    extern /* Subroutine */ int ehg182_(long *), ehg192_(double *, 
	    long *, long *, long *, long *, long *, double 
	    *, double *, long *);

    /* Parameter adjustments */
    --wv;
    --iv;
    --yy;

    /* Function Body */
    ++execnt;
    if (! (iv[28] != 172)) {
	ehg182_(&c__172);
    }
    if (! (iv[28] == 173)) {
	ehg182_(&c__173);
    }
    ehg192_(&yy[1], &iv[2], &iv[3], &iv[19], &iv[6], &iv[14], &wv[iv[13]], &
	    wv[iv[34]], &iv[iv[25]]);
    return 0;
} /* lowesr_ */

/* Subroutine */ int lowesw(double *res, long *n, double *rw, 
	long *pi)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    static double cmad;
    extern /* Subroutine */ int ehg106_(long *, long *, long *, 
	    long *, double *, long *, long *);
    static long i, i1;
    extern double d1mach_(long *);
    static long nh, identi;
    static double rsmall;
    extern long ifloor_(double *);

    /* Parameter adjustments */
    --pi;
    --rw;
    --res;

    /* Function Body */
/*     Identity -> identi */
    ++execnt;
/*     tranliterated from Devlin's ratfor */
/*     find median of absolute residuals */
    i__1 = *n;
    for (i1 = 1; i1 <= i__1; ++i1) {
	rw[i1] = (d__1 = res[i1], abs(d__1));
/* L3: */
    }
    i__1 = *n;
    for (identi = 1; identi <= i__1; ++identi) {
	pi[identi] = identi;
/* L4: */
    }
    d__1 = (double) (*n) / 2.;
    nh = ifloor_(&d__1) + 1;
/*     partial sort to find 6*mad */
    ehg106_(&c__1, n, &nh, &c__1, &rw[1], &pi[1], n);
    if (*n - nh + 1 < nh) {
	i__1 = nh - 1;
	i__2 = nh - 1;
	ehg106_(&c__1, &i__1, &i__2, &c__1, &rw[1], &pi[1], n);
	cmad = (rw[pi[nh]] + rw[pi[nh - 1]]) * 3;
    } else {
	cmad = rw[pi[nh]] * 6;
    }
    rsmall = d1mach_(&c__1);
    if (cmad < rsmall) {
	i__1 = *n;
	for (i1 = 1; i1 <= i__1; ++i1) {
	    rw[i1] = 1.;
/* L5: */
	}
    } else {
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    if (cmad * .999 < rw[i]) {
		rw[i] = 0.;
	    } else {
		if (cmad * .001 < rw[i]) {
/* Computing 2nd power */
		    d__2 = rw[i] / cmad;
/* Computing 2nd power */
		    d__1 = 1 - d__2 * d__2;
		    rw[i] = d__1 * d__1;
		} else {
		    rw[i] = 1.;
		}
	    }
/* L6: */
	}
    }
    return 0;
} /* lowesw_ */

/* Subroutine */ int lowesp(long *n, double *y, double *yhat, 
	double *pwgts, double *rwgts, long *pi, double *ytilde)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long i__1, i__2;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    extern /* Subroutine */ int ehg106_(long *, long *, long *, 
	    long *, double *, long *, long *);
    static double c;
    static long m;
    static double i1;
    static long i2, i3;
    static double i4;
    static long i5, identi;
    extern long ifloor_(double *);
    static double mad;

    /* Parameter adjustments */
    --ytilde;
    --pi;
    --rwgts;
    --pwgts;
    --yhat;
    --y;

    /* Function Body */
/*     Identity -> identi */
    ++execnt;
/*     median absolute deviation */
    i__1 = *n;
    for (i5 = 1; i5 <= i__1; ++i5) {
	ytilde[i5] = (d__1 = y[i5] - yhat[i5], abs(d__1)) * sqrt(pwgts[i5]);
/* L3: */
    }
    i__1 = *n;
    for (identi = 1; identi <= i__1; ++identi) {
	pi[identi] = identi;
/* L4: */
    }
    d__1 = (double) (*n) / 2.;
    m = ifloor_(&d__1) + 1;
    ehg106_(&c__1, n, &m, &c__1, &ytilde[1], &pi[1], n);
    if (*n - m + 1 < m) {
	i__1 = m - 1;
	i__2 = m - 1;
	ehg106_(&c__1, &i__1, &i__2, &c__1, &ytilde[1], &pi[1], n);
	mad = (ytilde[pi[m - 1]] + ytilde[pi[m]]) / 2;
    } else {
	mad = ytilde[pi[m]];
    }
/*     magic constant */
/* Computing 2nd power */
    d__1 = mad * 6;
    c = d__1 * d__1 / 5;
    i__1 = *n;
    for (i5 = 1; i5 <= i__1; ++i5) {
/* Computing 2nd power */
	d__1 = y[i5] - yhat[i5];
	ytilde[i5] = 1 - d__1 * d__1 * pwgts[i5] / c;
/* L5: */
    }
    i__1 = *n;
    for (i5 = 1; i5 <= i__1; ++i5) {
	ytilde[i5] *= sqrt(rwgts[i5]);
/* L6: */
    }
    if (*n <= 0) {
	i4 = 0.;
    } else {
	i3 = *n;
	i1 = ytilde[i3];
	for (i2 = i3 - 1; i2 >= 1; --i2) {
	    i1 = ytilde[i2] + i1;
/* L7: */
	}
	i4 = i1;
    }
    c = *n / i4;
/*     pseudovalues */
    i__1 = *n;
    for (i5 = 1; i5 <= i__1; ++i5) {
	ytilde[i5] = yhat[i5] + c * rwgts[i5] * (y[i5] - yhat[i5]);
/* L8: */
    }
    return 0;
} /* lowesp_ */

/* Subroutine */ int ehg124_(long *ll, long *uu, long *d, long *n,
	 long *nv, long *nc, long *ncmax, long *vc, double *x,
	 long *pi, long *a, double *xi, long *lo, long *hi, 
	long *c, double *v, long *vhit, long *nvmax, long *fc,
	 double *fd, long *dd)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long c_dim1, c_offset, v_dim1, v_offset, x_dim1, x_offset, i__1, i__2, 
	    i__3, i__4;
    double d__1;

    /* Builtin functions */
    double sqrt(double);
    long pow_ii(long *, long *);

    /* Local variables */
    static double diag[8];
    static long leaf;
    static double diam;
    extern /* Subroutine */ int ehg125_(long *, long *, double *, 
	    long *, long *, long *, long *, double *, long 
	    *, long *, long *, long *, long *), ehg106_(long *,
	     long *, long *, long *, double *, long *, 
	    long *), ehg129_(long *, long *, long *, double *,
	     long *, long *, double *);
    static long k, l, m, p, u;
    static double sigma[8];
    static long i1, i2, i3;
    static long i4, inorm2;
    extern long idamax_(long *, double *, long *);

    /* Parameter adjustments */
    --vhit;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --hi;
    --lo;
    --xi;
    --a;
    --pi;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    ++execnt;
    p = 1;
    l = *ll;
    u = *uu;
    lo[p] = l;
    hi[p] = u;
/*     top of while loop */
L3:
    if (! (p <= *nc)) {
	goto L4;
    }
    i__1 = *dd;
    for (i4 = 1; i4 <= i__1; ++i4) {
	diag[i4 - 1] = v[c[*vc + p * c_dim1] + i4 * v_dim1] - v[c[p * c_dim1 
		+ 1] + i4 * v_dim1];
/* L5: */
    }
    diam = 0.;
    i__1 = *dd;
    for (inorm2 = 1; inorm2 <= i__1; ++inorm2) {
/* Computing 2nd power */
	d__1 = diag[inorm2 - 1];
	diam += d__1 * d__1;
/* L6: */
    }
    diam = sqrt(diam);
    if (u - l + 1 <= *fc) {
	i1 = true;
    } else {
	i1 = diam <= *fd;
    }
    if (i1) {
	leaf = true;
    } else {
	if (*ncmax < *nc + 2) {
	    i2 = true;
	} else {
	    i2 = (double) (*nvmax) < *nv + (double) (*vc) / 2.;
	}
	leaf = i2;
    }
    if (! leaf) {
	ehg129_(&l, &u, dd, &x[x_offset], &pi[1], n, sigma);
	k = idamax_(dd, sigma, &c__1);
	m = (long) ((double) (l + u) / 2.);
	ehg106_(&l, &u, &m, &c__1, &x[k * x_dim1 + 1], &pi[1], n);
/*           all ties go with hi son */
/*           top of while loop */
L7:
	if (1 < m) {
	    i3 = x[pi[m - 1] + k * x_dim1] == x[pi[m] + k * x_dim1];
	} else {
	    i3 = false;
	}
	if (! i3) {
	    goto L8;
	}
	--m;
	goto L7;
/*           bottom of while loop */
L8:
	if (v[c[p * c_dim1 + 1] + k * v_dim1] == x[pi[m] + k * x_dim1]) {
	    leaf = true;
	} else {
	    leaf = v[c[*vc + p * c_dim1] + k * v_dim1] == x[pi[m] + k * 
		    x_dim1];
	}
    }
    if (leaf) {
	a[p] = 0;
    } else {
	a[p] = k;
	xi[p] = x[pi[m] + k * x_dim1];
/*           left son */
	++(*nc);
	lo[p] = *nc;
	lo[*nc] = l;
	hi[*nc] = m;
/*           right son */
	++(*nc);
	hi[p] = *nc;
	lo[*nc] = m + 1;
	hi[*nc] = u;
	i__2 = k - 1;
	i__1 = pow_ii(&c__2, &i__2);
	i__4 = *d - k;
	i__3 = pow_ii(&c__2, &i__4);
	ehg125_(&p, nv, &v[v_offset], &vhit[1], nvmax, d, &k, &xi[p], &i__1, &
		i__3, &c[p * c_dim1 + 1], &c[lo[p] * c_dim1 + 1], &c[hi[p] * 
		c_dim1 + 1]);
    }
    ++p;
    l = lo[p];
    u = hi[p];
    goto L3;
/*     bottom of while loop */
L4:
    return 0;
} /* ehg124_ */

/* Subroutine */ int ehg129_(long *l, long *u, long *d, double *
	x, long *pi, long *n, double *sigma)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long x_dim1, x_offset, i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    static double beta;
    static long i, k;
    static double t, alpha;
    extern double d1mach_(long *);
    static double machin;

    /* Parameter adjustments */
    --sigma;
    --pi;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
/*     MachInf -> machin */
    ++execnt;
    if (execnt == 1) {
	machin = d1mach_(&c__2);
    }
    i__1 = *d;
    for (k = 1; k <= i__1; ++k) {
	alpha = machin;
	beta = -machin;
	i__2 = *u;
	for (i = *l; i <= i__2; ++i) {
	    t = x[pi[i] + k * x_dim1];
/* Computing MIN */
	    d__1 = alpha, d__2 = x[pi[i] + k * x_dim1];
	    alpha = min(d__1,d__2);
	    beta = max(beta,t);
/* L4: */
	}
	sigma[k] = beta - alpha;
/* L3: */
    }
    return 0;
} /* ehg129_ */

/* Subroutine */ int ehg137_(double *z, long *kappa, long *leaf, 
	long *nleaf, long *d, long *nv, long *nvmax, long *
	ncmax, long *vc, long *a, double *xi, long *lo, long *
	hi, long *c, double *v)
{
    /* Initialized data */

    static long execnt = 0;

    /* System generated locals */
    long i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int ehg182_(long *);
    static long p, pstack[20], stackt;

    /* Parameter adjustments */
    --hi;
    --lo;
    --xi;
    --a;
    --leaf;
    --z;

    /* Function Body */
/*     stacktop -> stackt */
    ++execnt;
/*     find leaf cells affected by $z$ */
    stackt = 0;
    p = 1;
    *nleaf = 0;
/*     top of while loop */
L3:
    if (! (0 < p)) {
	goto L4;
    }
    if (a[p] == 0) {
/*           leaf */
	++(*nleaf);
	leaf[*nleaf] = p;
/*           Pop */
	if (stackt >= 1) {
	    p = pstack[stackt - 1];
	} else {
	    p = 0;
	}
/* Computing MAX */
	i__1 = 0, i__2 = stackt - 1;
	stackt = max(i__1,i__2);
    } else {
	if (z[a[p]] == xi[p]) {
/*              Push */
	    ++stackt;
	    if (! (stackt <= 20)) {
		ehg182_(&c__187);
	    }
	    pstack[stackt - 1] = hi[p];
	    p = lo[p];
	} else {
	    if (z[a[p]] < xi[p]) {
		p = lo[p];
	    } else {
		p = hi[p];
	    }
	}
    }
    goto L3;
/*     bottom of while loop */
L4:
    if (! (*nleaf <= 256)) {
	ehg182_(&c__185);
    }
    return 0;
} /* ehg137_ */


// Fortran Functions

double gamma(double xx)
{
   int j;
   double x,y,tmp,ser;
   double cof[6] = {
      76.18009172947146,    -86.50532032941677,
      24.01409824083091,    -1.231739572450155,
      0.1208650973866179e-2,-0.5395239384953e-5
   };

   y = x = xx;
   tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
   ser = 1.000000000190015;
   for (j=0;j<=5;j++)
      ser += (cof[j] / ++y);
   return log(2.5066282746310005 * ser / x) - tmp;
}

/*  If your compiler is so ancient it doesn't recognize void, say
#define void
*/

void anova(struct loess_struct *one, struct loess_struct *two, struct anova_struct *out)
{
	double	one_d1, one_d2, one_s, two_d1, two_d2, two_s,
	        rssdiff, d1diff, tmp;
	int     max_enp;

	one_d1 = one->out.one_delta;
	one_d2 = one->out.two_delta;
	one_s = one->out.s;
	two_d1 = two->out.one_delta;
	two_d2 = two->out.two_delta;
	two_s = two->out.s;

        rssdiff = fabs(one_s * one_s * one_d1 - two_s * two_s * two_d1);
        d1diff = fabs(one_d1 - two_d1);
        out->dfn = d1diff * d1diff / fabs(one_d2 - two_d2);
	max_enp = (one->out.enp > two->out.enp);
	tmp = max_enp ? one_d1 : two_d1;
        out->dfd = tmp * tmp / (max_enp ? one_d2 : two_d2);
	tmp = max_enp ? one_s : two_s;
        out->F_value = (rssdiff / d1diff) / (tmp * tmp);
        out->Pr_F = 1 - pf(out->F_value, out->dfn, out->dfd);
}

void pointwise(struct  pred_struct *pre, long m, double coverage, struct  ci_struct *ci)
{	
	double	t_dist, limit, fit;
	int	i;	

        ci->fit = (double *) malloc(m * sizeof(double));
        ci->upper = (double *) malloc(m * sizeof(double));
	ci->lower = (double *) malloc(m * sizeof(double));

	t_dist = qt(1 - (1 - coverage)/2, pre->df);
	for(i = 0; i < m; i++) {
		limit = pre->se_fit[i] * t_dist;
		ci->fit[i] = fit = pre->fit[i];
		ci->upper[i] = fit + limit;
		ci->lower[i] = fit - limit;
	}	
}

void pw_free_mem(struct ci_struct *ci)
{
    free(ci->fit);
    free(ci->upper);
	free(ci->lower);
}

double pf(double q, double df1, double df2)
{
	return(ibeta(q*df1/(df2+q*df1), df1/2, df2/2));
}

double qt(double p, double df)
{
        double        t;

	t = invibeta(fabs(2*p-1), 0.5, df/2);
        return((p>0.5?1:-1) * sqrt(t*df/(1-t)));
}

/**********************************************************************/
 /*
 * Incomplete beta function.
 * Reference:  Abramowitz and Stegun, 26.5.8.
 * Assumptions: 0 <= x <= 1; a,b > 0.
 */
#define DOUBLE_EPS      2.2204460492503131E-16
#define IBETA_LARGE     1.0e30
#define IBETA_SMALL     1.0e-30

double ibeta(double x, double a, double b)
{
        int flipped = 0, i, k, count;
        double I, temp, pn[6], ak, bk, next, prev, factor, val;

        if (x <= 0)
                return(0);
        if (x >= 1)
                return(1);

        /* use ibeta(x,a,b) = 1-ibeta(1-x,b,a) */
        if ((a+b+1)*x > (a+1)) {
                flipped = 1;
                temp = a;
                a = b;
                b = temp;
                x = 1 - x;
        }

        pn[0] = 0.0;
        pn[2] = pn[3] = pn[1] = 1.0;
        count = 1;
        val = x/(1.0-x);
        bk = 1.0;
        next = 1.0;
        do {
                count++;
                k = count/2;
                prev = next;
                if (count%2 == 0)
                        ak = -((a+k-1.0)*(b-k)*val)/
                                ((a+2.0*k-2.0)*(a+2.0*k-1.0));
                else
                        ak = ((a+b+k-1.0)*k*val)/
                                ((a+2.0*k)*(a+2.0*k-1.0));
                pn[4] = bk*pn[2] + ak*pn[0];
                pn[5] = bk*pn[3] + ak*pn[1];
                next = pn[4] / pn[5];
                for (i=0; i<=3; i++)
                        pn[i] = pn[i+2];
                if (fabs(pn[4]) >= IBETA_LARGE)
                        for (i=0; i<=3; i++)
                                pn[i] /= IBETA_LARGE;
                if (fabs(pn[4]) <= IBETA_SMALL)
                        for (i=0; i<=3; i++)
                                pn[i] /= IBETA_SMALL;
        } while (fabs(next-prev) > DOUBLE_EPS*prev);
        factor = a*log(x) + (b-1)*log(1-x);
        factor -= gamma(a+1) + gamma(b) - gamma(a+b);
        I = exp(factor) * next;
        return(flipped ? 1-I : I);
}

/*
 * Rational approximation to inverse Gaussian distribution.
 * Absolute error is bounded by 4.5e-4.
 * Reference: Abramowitz and Stegun, page 933.
 * Assumption: 0 < p < 1.
 */

static double num[] = {
        2.515517,
        0.802853,
        0.010328
};

static double den[] = {
        1.000000,
        1.432788,
        0.189269,
        0.001308
};

double invigauss_quick(double p)
{
        int lower;
        double t, n, d, q;

        if(p == 0.5)
                return(0);
        lower = p < 0.5;
        p = lower ? p : 1 - p;
        t = sqrt(-2 * log(p));
        n = (num[2]*t + num[1])*t + num[0];
        d = ((den[3]*t + den[2])*t + den[1])*t + den[0];
        q = lower ? n/d - t : t - n/d;
        return(q);
}

/*
 * Inverse incomplete beta function.
 * Assumption: 0 <= p <= 1, a,b > 0.
 */

double invibeta(double p, double a, double b)
{
        int i;
        double ql, qr, qm, qdiff;
        double pl, pr, pm, pdiff;

/*        MEANINGFUL(qm);*/
	qm = 0;
        if(p == 0 || p == 1)
                return(p);

        /* initialize [ql,qr] containing the root */
        ql = qr = invibeta_quick(p, a, b);
        pl = pr = ibeta(ql, a, b);
        if(pl == p)
                return(ql);
        if(pl < p)
                while(1) {
                        qr += 0.05;
                        if(qr >= 1) {
                                pr = qr = 1;
                                break;
                        }
                        pr = ibeta(qr, a, b);
                        if(pr == p)
                                return(pr);
                        if(pr > p)
                                break;
                }
        else
                while(1) {
                        ql -= 0.05;
                        if(ql <= 0) {
                                pl = ql = 0;
                                break;
                        }
                        pl = ibeta(ql, a, b);
                        if(pl == p)
                                return(pl);
                        if(pl < p)
                                break;
                }

        /* a few steps of bisection */
        for(i = 0; i < 5; i++) {
                qm = (ql + qr) / 2;
                pm = ibeta(qm, a, b);
                qdiff = qr - ql;
                pdiff = pm - p;
                if(fabs(qdiff) < DOUBLE_EPS*qm || fabs(pdiff) < DOUBLE_EPS)
                        return(qm);
                if(pdiff < 0) {
                        ql = qm;
                        pl = pm;
                } else {
                        qr = qm;
                        pr = pm;
                }
        }

        /* a few steps of secant */
        for(i = 0; i < 40; i++) {
                qm = ql + (p-pl)*(qr-ql)/(pr-pl);
                pm = ibeta(qm, a, b);
                qdiff = qr - ql;
                pdiff = pm - p;
                if(fabs(qdiff) < 2*DOUBLE_EPS*qm || fabs(pdiff) < 2*DOUBLE_EPS)
                        return(qm);
                if(pdiff < 0) {
                        ql = qm;
                        pl = pm;
                } else {
                        qr = qm;
                        pr = pm;
                }
        }

        /* no convergence */
        return(qm);
}

/*
 * Quick approximation to inverse incomplete beta function,
 * by matching first two moments with the Gaussian distribution.
 * Assumption: 0 < p < 1, a,b > 0.
 */

double invibeta_quick(double p, double a, double b)
{
        double x, m, s;

        x = a + b;
        m = a / x;
        s = sqrt((a*b) / (x*x*(x+1)));
        return(fmax(0.0, fmin(1.0, invigauss_quick(p)*s + m)));
}
 
static double fmin(double a, double b)
{
        return(a < b ? a : b);
}

static double fmax(double a, double b)
{
        return(a > b ? a : b);
}

void Recover(char* a, int *b)
{
        printf(a);
        exit(1);
}

void Warning(char *a, int *b)
{
        printf(a);
}

/*  d1mach may be replaced by Fortran code:
    mail netlib@netlib.bell-labs.com
    send d1mach from core.
*/

double d1mach_ ( long *i)
{
	switch(*i){
	case 1: return DBL_MIN;
	case 2: return DBL_MAX;
	case 3: return DBL_EPSILON/FLT_RADIX;
	case 4: return DBL_EPSILON;
	case 5: return log10(FLT_RADIX);
        default: Recover("Invalid argument to d1mach()", 0L);
        }

	return 0.0;
}

double ddot_(long *n, double *dx, long *incx, double *dy, long *incy)
{


    /* System generated locals */
    long i__1;
    double ret_val;

    /* Local variables */
    static long i, m;
    static double dtemp;
    static long ix, iy, mp1;


/*     forms the dot product of two vectors.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp += DX(ix) * DY(iy);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	dtemp += DX(i) * DY(i);
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 5) {
	dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) * 
		DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

/* Subroutine */ int dqrsl_(
double *x, 
long *ldx, 
long *n, 
long *k, 
double *qraux, 
double *y, 
double *qy, 
double *qty, 
double *b, 
double *rsd, 
double *xb, 
long *job, 
long *info)
{
    /* System generated locals */
    long x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static double temp;
    static long cqty;
    static long i, j;
    static double t;
    static long cb;
    static long jj;
    static long cr;
    static long ju, kp1;
    static long cxb, cqy;

/* ***BEGIN PROLOGUE  DQRSL */
/* ***DATE WRITTEN   780814   (YYMMDD) */
/* ***REVISION DATE  820801   (YYMMDD) */
/* ***CATEGORY NO.  D9,D2A1 */
/* ***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX, */
/*             ORTHOGONAL TRIANGULAR,SOLVE */
/* ***AUTHOR  STEWART, G. W., (U. OF MARYLAND) */
/* ***PURPOSE  Applies the output of DQRDC to compute coordinate */
/*            transformations, projections, and least squares solutions. 
*/
/* ***DESCRIPTION */

/*     DQRSL applies the output of DQRDC to compute coordinate */
/*     transformations, projections, and least squares solutions. */
/*     For K .LE. MIN(N,P), let XK be the matrix */

/*            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K))) */

/*     formed from columnns JPVT(1), ... ,JPVT(K) of the original */
/*     N X P matrix X that was input to DQRDC (if no pivoting was */
/*     done, XK consists of the first K columns of X in their */
/*     original order).  DQRDC produces a factored orthogonal matrix Q */
/*     and an upper triangular matrix R such that */

/*              XK = Q * (R) */
/*                       (0) */

/*     This information is contained in coded form in the arrays */
/*     X and QRAUX. */

/*     On Entry */

/*        X      DOUBLE PRECISION(LDX,P). */
/*               X contains the output of DQRDC. */

/*        LDX    long. */
/*               LDX is the leading dimension of the array X. */

/*        N      long. */
/*               N is the number of rows of the matrix XK.  It must */
/*               have the same value as N in DQRDC. */

/*        K      long. */
/*               K is the number of columns of the matrix XK.  K */
/*               must not be greater than MIN(N,P), where P is the */
/*               same as in the calling sequence to DQRDC. */

/*        QRAUX  DOUBLE PRECISION(P). */
/*               QRAUX contains the auxiliary output from DQRDC. */

/*        Y      DOUBLE PRECISION(N) */
/*               Y contains an N-vector that is to be manipulated */
/*               by DQRSL. */

/*        JOB    long. */
/*               JOB specifies what is to be computed.  JOB has */
/*               the decimal expansion ABCDE, with the following */
/*               meaning. */

/*                    If A .NE. 0, compute QY. */
/*                    If B,C,D, or E .NE. 0, compute QTY. */
/*                    If C .NE. 0, compute B. */
/*                    If D .NE. 0, compute RSD. */
/*                    If E .NE. 0, compute XB. */

/*               Note that a request to compute B, RSD, or XB */
/*               automatically triggers the computation of QTY, for */
/*               which an array must be provided in the calling */
/*               sequence. */

/*     On Return */

/*        QY     DOUBLE PRECISION(N). */
/*               QY contains Q*Y, if its computation has been */
/*               requested. */

/*        QTY    DOUBLE PRECISION(N). */
/*               QTY contains TRANS(Q)*Y, if its computation has */
/*               been requested.  Here TRANS(Q) is the */
/*               transpose of the matrix Q. */

/*        B      DOUBLE PRECISION(K) */
/*               B contains the solution of the least squares problem */

/*                    minimize norm2(Y - XK*B), */

/*               if its computation has been requested.  (Note that */
/*               if pivoting was requested in DQRDC, the J-th */
/*               component of B will be associated with column JPVT(J) */
/*               of the original matrix X that was input into DQRDC.) */

/*        RSD    DOUBLE PRECISION(N). */
/*               RSD contains the least squares residual Y - XK*B, */
/*               if its computation has been requested.  RSD is */
/*               also the orthogonal projection of Y onto the */
/*               orthogonal complement of the column space of XK. */

/*        XB     DOUBLE PRECISION(N). */
/*               XB contains the least squares approximation XK*B, */
/*               if its computation has been requested.  XB is also */
/*               the orthogonal projection of Y onto the column space */
/*               of X. */

/*        INFO   long. */
/*               INFO is zero unless the computation of B has */
/*               been requested and R is exactly singular.  In */
/*               this case, INFO is the index of the first zero */
/*               diagonal element of R and B is left unaltered. */

/*     The parameters QY, QTY, B, RSD, and XB are not referenced */
/*     if their computation is not requested and in this case */
/*     can be replaced by dummy variables in the calling program. */
/*     To save storage, the user may in some cases use the same */
/*     array for different parameters in the calling sequence.  A */
/*     frequently occuring example is when one wishes to compute */
/*     any of B, RSD, or XB and does not need Y or QTY.  In this */
/*     case one may identify Y, QTY, and one of B, RSD, or XB, while */
/*     providing separate arrays for anything else that is to be */
/*     computed.  Thus the calling sequence */

/*          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO) */

/*     will result in the computation of B and RSD, with RSD */
/*     overwriting Y.  More generally, each item in the following */
/*     list contains groups of permissible identifications for */
/*     a single calling sequence. */

/*          1. (Y,QTY,B) (RSD) (XB) (QY) */

/*          2. (Y,QTY,RSD) (B) (XB) (QY) */

/*          3. (Y,QTY,XB) (B) (RSD) (QY) */

/*          4. (Y,QY) (QTY,B) (RSD) (XB) */

/*          5. (Y,QY) (QTY,RSD) (B) (XB) */

/*          6. (Y,QY) (QTY,XB) (B) (RSD) */

/*     In any group the value returned in the array allocated to */
/*     the group corresponds to the last member of the group. */

/*     LINPACK.  This version dated 08/14/78 . */
/*     G. W. Stewart, University of Maryland, Argonne National Lab. */

/*     DQRSL uses the following functions and subprograms. */

/*     BLAS DAXPY,DCOPY,DDOT */
/*     Fortran DABS,MIN0,MOD */
/* ***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W., */
/*                 *LINPACK USERS  GUIDE*, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY,DCOPY,DDOT */
/* ***END PROLOGUE  DQRSL */


/*     SET INFO long. */

/* ***FIRST EXECUTABLE STATEMENT  DQRSL */
    /* Parameter adjustments */
    --xb;
    --rsd;
    --b;
    --qty;
    --qy;
    --y;
    --qraux;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    *info = 0;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    cqy = *job / 10000 != 0;
    cqty = *job % 10000 != 0;
    cb = *job % 1000 / 100 != 0;
    cr = *job % 100 / 10 != 0;
    cxb = *job % 10 != 0;
/* Computing MIN */
    i__1 = *k, i__2 = *n - 1;
    ju = min(i__1,i__2);

/*     SPECIAL ACTION WHEN N=1. */

    if (ju != 0) {
	goto L40;
    }
    if (cqy) {
	qy[1] = y[1];
    }
    if (cqty) {
	qty[1] = y[1];
    }
    if (cxb) {
	xb[1] = y[1];
    }
    if (! cb) {
	goto L30;
    }
    if (x[x_dim1 + 1] != 0.) {
	goto L10;
    }
    *info = 1;
    goto L20;
L10:
    b[1] = y[1] / x[x_dim1 + 1];
L20:
L30:
    if (cr) {
	rsd[1] = 0.;
    }
    goto L250;
L40:

/*        SET UP TO COMPUTE QY OR QTY. */

    if (cqy) {
	dcopy_(n, &y[1], &c__1, &qy[1], &c__1);
    }
    if (cqty) {
	dcopy_(n, &y[1], &c__1, &qty[1], &c__1);
    }
    if (! cqy) {
	goto L70;
    }

/*           COMPUTE QY. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	if (qraux[j] == 0.) {
	    goto L50;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &qy[j], &c__1) / x[j + j 
		* x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qy[j], &c__1);
	x[j + j * x_dim1] = temp;
L50:
/* L60: */
	;
    }
L70:
    if (! cqty) {
	goto L100;
    }

/*           COMPUTE TRANS(Q)*Y. */

    i__1 = ju;
    for (j = 1; j <= i__1; ++j) {
	if (qraux[j] == 0.) {
	    goto L80;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &qty[j], &c__1) / x[j + 
		j * x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
	x[j + j * x_dim1] = temp;
L80:
/* L90: */
	;
    }
L100:

/*        SET UP TO COMPUTE B, RSD, OR XB. */

    if (cb) {
	dcopy_(k, &qty[1], &c__1, &b[1], &c__1);
    }
    kp1 = *k + 1;
    if (cxb) {
	dcopy_(k, &qty[1], &c__1, &xb[1], &c__1);
    }
    if (cr && *k < *n) {
	i__1 = *n - *k;
	dcopy_(&i__1, &qty[kp1], &c__1, &rsd[kp1], &c__1);
    }
    if (! cxb || kp1 > *n) {
	goto L120;
    }
    i__1 = *n;
    for (i = kp1; i <= i__1; ++i) {
	xb[i] = 0.;
/* L110: */
    }
L120:
    if (! cr) {
	goto L140;
    }
    i__1 = *k;
    for (i = 1; i <= i__1; ++i) {
	rsd[i] = 0.;
/* L130: */
    }
L140:
    if (! cb) {
	goto L190;
    }

/*           COMPUTE B. */

    i__1 = *k;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *k - jj + 1;
	if (x[j + j * x_dim1] != 0.) {
	    goto L150;
	}
	*info = j;
/*           ......EXIT */
	goto L180;
L150:
	b[j] /= x[j + j * x_dim1];
	if (j == 1) {
	    goto L160;
	}
	t = -b[j];
	i__2 = j - 1;
	daxpy_(&i__2, &t, &x[j * x_dim1 + 1], &c__1, &b[1], &c__1);
L160:
/* L170: */
	;
    }
L180:
L190:
    if (! cr && ! cxb) {
	goto L240;
    }

/*           COMPUTE RSD OR XB AS REQUIRED. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	if (qraux[j] == 0.) {
	    goto L220;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	if (! cr) {
	    goto L200;
	}
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1) / x[j + 
		j * x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
L200:
	if (! cxb) {
	    goto L210;
	}
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &xb[j], &c__1) / x[j + j 
		* x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
L210:
	x[j + j * x_dim1] = temp;
L220:
/* L230: */
	;
    }
L240:
L250:
    return 0;
} /* dqrsl_ */

int dcopy_(long *n, double *dx, long *incx, 
	double *dy, long *incy)
{


    /* System generated locals */
    long i__1;

    /* Local variables */
    static long i, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	DY(iy) = DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	DY(i) = DX(i);
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 7) {
	DY(i) = DX(i);
	DY(i + 1) = DX(i + 1);
	DY(i + 2) = DX(i + 2);
	DY(i + 3) = DX(i + 3);
	DY(i + 4) = DX(i + 4);
	DY(i + 5) = DX(i + 5);
	DY(i + 6) = DX(i + 6);
/* L50: */
    }
    return 0;
} /* dcopy_ */

int daxpy_(long *n, double *da, double *dx, 
	long *incx, double *dy, long *incy)
{


    /* System generated locals */
    long i__1;

    /* Local variables */
    static long i, m, ix, iy, mp1;


/*     constant times a vector plus a vector.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	DY(iy) += *da * DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	DY(i) += *da * DX(i);
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
	DY(i) += *da * DX(i);
	DY(i + 1) += *da * DX(i + 1);
	DY(i + 2) += *da * DX(i + 2);
	DY(i + 3) += *da * DX(i + 3);
/* L50: */
    }
    return 0;
} /* daxpy_ */

int dsvdc_(double *x, long *ldx, long *n, long *p, double *s, double *e, double *u, long *ldu, double *v, long *ldv, double *work, long *job, long *info)
{
    /* System generated locals */
    long x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;
    double d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Local variables */
    static long kase;
    static long jobu, iter;
    static double test;
    static long nctp1;
    static double b, c;
    static long nrtp1;
    static double f, g;
    static long i, j, k, l, m;
    static double t, scale;
    static double shift;
    static long maxit;
    static long wantu, wantv;
    static double t1, ztest, el;
    static long kk;
    static double cs;
    static long ll, mm, ls;
    static double sl;
    static long lu;
    static double sm, sn;
    static long lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    static double emm1, smm1;

/* ***BEGIN PROLOGUE  DSVDC */
/* ***DATE WRITTEN   790319   (YYMMDD) */
/* ***REVISION DATE  820801   (YYMMDD) */
/* ***CATEGORY NO.  D6 */
/* ***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX, */
/*             SINGULAR VALUE DECOMPOSITION */
/* ***AUTHOR  STEWART, G. W., (U. OF MARYLAND) */
/* ***PURPOSE  Perform the singular value decomposition of a d.p. NXP */
/*            matrix. */
/* ***DESCRIPTION */

/*     DSVDC is a subroutine to reduce a double precision NxP matrix X */
/*     by orthogonal transformations U and V to diagonal form.  The */
/*     diagonal elements S(I) are the singular values of X.  The */
/*     columns of U are the corresponding left singular vectors, */
/*     and the columns of V the right singular vectors. */

/*     On Entry */

/*         X         DOUBLE PRECISION(LDX,P), where LDX .GE. N. */
/*                   X contains the matrix whose singular value */
/*                   decomposition is to be computed.  X is */
/*                   destroyed by DSVDC. */

/*         LDX       long. */
/*                   LDX is the leading dimension of the array X. */

/*         N         long. */
/*                   N is the number of columns of the matrix X. */

/*         P         long. */
/*                   P is the number of rows of the matrix X. */

/*         LDU       long. */
/*                   LDU is the leading dimension of the array U. */
/*                   (See below). */

/*         LDV       long. */
/*                   LDV is the leading dimension of the array V. */
/*                   (See below). */

/*         WORK      DOUBLE PRECISION(N). */
/*                   WORK is a scratch array. */

/*         JOB       long. */
/*                   JOB controls the computation of the singular */
/*                   vectors.  It has the decimal expansion AB */
/*                   with the following meaning */

/*                        A .EQ. 0    do not compute the left singular */
/*                                  vectors. */
/*                        A .EQ. 1    return the N left singular vectors 
*/
/*                                  in U. */
/*                        A .GE. 2    return the first MIN(N,P) singular 
*/
/*                                  vectors in U. */
/*                        B .EQ. 0    do not compute the right singular */
/*                                  vectors. */
/*                        B .EQ. 1    return the right singular vectors */
/*                                  in V. */

/*     On Return */

/*         S         DOUBLE PRECISION(MM), where MM=MIN(N+1,P). */
/*                   The first MIN(N,P) entries of S contain the */
/*                   singular values of X arranged in descending */
/*                   order of magnitude. */

/*         E         DOUBLE PRECISION(P). */
/*                   E ordinarily contains zeros.  However see the */
/*                   discussion of INFO for exceptions. */

/*         U         DOUBLE PRECISION(LDU,K), where LDU .GE. N. */
/*                   If JOBA .EQ. 1, then K .EQ. N. */
/*                   If JOBA .GE. 2, then K .EQ. MIN(N,P). */
/*                   U contains the matrix of right singular vectors. */
/*                   U is not referenced if JOBA .EQ. 0.  If N .LE. P */
/*                   or if JOBA .EQ. 2, then U may be identified with X */
/*                   in the subroutine call. */

/*         V         DOUBLE PRECISION(LDV,P), where LDV .GE. P. */
/*                   V contains the matrix of right singular vectors. */
/*                   V is not referenced if JOB .EQ. 0.  If P .LE. N, */
/*                   then V may be identified with X in the */
/*                   subroutine call. */

/*         INFO      long. */
/*                   The singular values (and their corresponding */
/*                   singular vectors) S(INFO+1),S(INFO+2),...,S(M) */
/*                   are correct (here M=MIN(N,P)).  Thus if */
/*                   INFO .EQ. 0, all the singular values and their */
/*                   vectors are correct.  In any event, the matrix */
/*                   B = TRANS(U)*X*V is the bidiagonal matrix */
/*                   with the elements of S on its diagonal and the */
/*                   elements of E on its super-diagonal (TRANS(U) */
/*                   is the transpose of U).  Thus the singular */
/*                   values of X and B are the same. */

/*     LINPACK.  This version dated 03/19/79 . */
/*     G. W. Stewart, University of Maryland, Argonne National Lab. */

/*     DSVDC uses the following functions and subprograms. */

/*     External DROT */
/*     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG */
/*     Fortran DABS,DMAX1,MAX0,MIN0,MOD,DSQRT */
/* ***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W., */
/*                 *LINPACK USERS  GUIDE*, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY,DDOT,DNRM2,DROT,DROTG,DSCAL,DSWAP */
/* ***END PROLOGUE  DSVDC */



/*     SET THE MAXIMUM NUMBER OF ITERATIONS. */

/* ***FIRST EXECUTABLE STATEMENT  DSVDC */
    /* Parameter adjustments */
    --work;
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    u_dim1 = *ldu;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --e;
    --s;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    maxit = 30;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    wantu = false;
    wantv = false;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1) {
	ncu = min(*n,*p);
    }
    if (jobu != 0) {
	wantu = true;
    }
    if (*job % 10 != 0) {
	wantv = true;
    }

/*     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS */
/*     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E. */

    *info = 0;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    if (lu < 1) {
	goto L170;
    }
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
	lp1 = l + 1;
	if (l > nct) {
	    goto L20;
	}

/*           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND */
/*           PLACE THE L-TH DIAGONAL IN S(L). */

	i__2 = *n - l + 1;
	s[l] = dnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (s[l] == 0.) {
	    goto L10;
	}
	if (x[l + l * x_dim1] != 0.) {
	    s[l] = d_sign(&s[l], &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / s[l];
	dscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;
L10:
	s[l] = -s[l];
L20:
	if (*p < lp1) {
	    goto L50;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    if (l > nct) {
		goto L30;
	    }
	    if (s[l] == 0.) {
		goto L30;
	    }

/*              APPLY THE TRANSFORMATION. */

	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
L30:

/*           PLACE THE L-TH ROW OF X INTO  E FOR THE */
/*           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION. */

	    e[j] = x[l + j * x_dim1];
/* L40: */
	}
L50:
	if (! wantu || l > nct) {
	    goto L70;
	}

/*           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK */
/*           MULTIPLICATION. */

	i__2 = *n;
	for (i = l; i <= i__2; ++i) {
	    u[i + l * u_dim1] = x[i + l * x_dim1];
/* L60: */
	}
L70:
	if (l > nrt) {
	    goto L150;
	}

/*           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE */
/*           L-TH SUPER-DIAGONAL IN E(L). */

	i__2 = *p - l;
	e[l] = dnrm2_(&i__2, &e[lp1], &c__1);
	if (e[l] == 0.) {
	    goto L80;
	}
	if (e[lp1] != 0.) {
	    e[l] = d_sign(&e[l], &e[lp1]);
	}
	i__2 = *p - l;
	d__1 = 1. / e[l];
	dscal_(&i__2, &d__1, &e[lp1], &c__1);
	e[lp1] += 1.;
L80:
	e[l] = -e[l];
	if (lp1 > *n || e[l] == 0.) {
	    goto L120;
	}

/*              APPLY THE TRANSFORMATION. */

	i__2 = *n;
	for (i = lp1; i <= i__2; ++i) {
	    work[i] = 0.;
/* L90: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    daxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &
		    c__1);
/* L100: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    d__1 = -e[j] / e[lp1];
	    daxpy_(&i__3, &d__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &
		    c__1);
/* L110: */
	}
L120:
	if (! wantv) {
	    goto L140;
	}

/*              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT */
/*              BACK MULTIPLICATION. */

	i__2 = *p;
	for (i = lp1; i <= i__2; ++i) {
	    v[i + l * v_dim1] = e[i];
/* L130: */
	}
L140:
L150:
/* L160: */
	;
    }
L170:

/*     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M. */

/* Computing MIN */
    i__1 = *p, i__2 = *n + 1;
    m = min(i__1,i__2);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p) {
	s[nctp1] = x[nctp1 + nctp1 * x_dim1];
    }
    if (*n < m) {
	s[m] = 0.;
    }
    if (nrtp1 < m) {
	e[nrtp1] = x[nrtp1 + m * x_dim1];
    }
    e[m] = 0.;

/*     IF REQUIRED, GENERATE U. */

    if (! wantu) {
	goto L300;
    }
    if (ncu < nctp1) {
	goto L200;
    }
    i__1 = ncu;
    for (j = nctp1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    u[i + j * u_dim1] = 0.;
/* L180: */
	}
	u[j + j * u_dim1] = 1.;
/* L190: */
    }
L200:
    if (nct < 1) {
	goto L290;
    }
    i__1 = nct;
    for (ll = 1; ll <= i__1; ++ll) {
	l = nct - ll + 1;
	if (s[l] == 0.) {
	    goto L250;
	}
	lp1 = l + 1;
	if (ncu < lp1) {
	    goto L220;
	}
	i__2 = ncu;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1) / u[l + l * u_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1);
/* L210: */
	}
L220:
	i__2 = *n - l + 1;
	dscal_(&i__2, &c_b44, &u[l + l * u_dim1], &c__1);
	u[l + l * u_dim1] += 1.;
	lm1 = l - 1;
	if (lm1 < 1) {
	    goto L240;
	}
	i__2 = lm1;
	for (i = 1; i <= i__2; ++i) {
	    u[i + l * u_dim1] = 0.;
/* L230: */
	}
L240:
	goto L270;
L250:
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    u[i + l * u_dim1] = 0.;
/* L260: */
	}
	u[l + l * u_dim1] = 1.;
L270:
/* L280: */
	;
    }
L290:
L300:

/*     IF IT IS REQUIRED, GENERATE V. */

    if (! wantv) {
	goto L350;
    }
    i__1 = *p;
    for (ll = 1; ll <= i__1; ++ll) {
	l = *p - ll + 1;
	lp1 = l + 1;
	if (l > nrt) {
	    goto L320;
	}
	if (e[l] == 0.) {
	    goto L320;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *p - l;
	    t = -ddot_(&i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1) / v[lp1 + l * v_dim1];
	    i__3 = *p - l;
	    daxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
/* L310: */
	}
L320:
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    v[i + l * v_dim1] = 0.;
/* L330: */
	}
	v[l + l * v_dim1] = 1.;
/* L340: */
    }
L350:

/*     MAIN ITERATION LOOP FOR THE SINGULAR VALUES. */

    mm = m;
    iter = 0;
L360:

/*        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND. */

/*     ...EXIT */
    if (m == 0) {
	goto L620;
    }

/*        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET */
/*        long AND RETURN. */

    if (iter < maxit) {
	goto L370;
    }
    *info = m;
/*     ......EXIT */
    goto L620;
L370:

/*        THIS SECTION OF THE PROGRAM INSPECTS FOR */
/*        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON */
/*        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS. */

/*           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M */
/*           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M */
/*           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND */
/*                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP). */
/*           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE). */

    i__1 = m;
    for (ll = 1; ll <= i__1; ++ll) {
	l = m - ll;
/*        ...EXIT */
	if (l == 0) {
	    goto L400;
	}
	test = (d__1 = s[l], abs(d__1)) + (d__2 = s[l + 1], abs(d__2));
	ztest = test + (d__1 = e[l], abs(d__1));
	if (ztest != test) {
	    goto L380;
	}
	e[l] = 0.;
/*        ......EXIT */
	goto L400;
L380:
/* L390: */
	;
    }
L400:
    if (l != m - 1) {
	goto L410;
    }
    kase = 4;
    goto L480;
L410:
    lp1 = l + 1;
    mp1 = m + 1;
    i__1 = mp1;
    for (lls = lp1; lls <= i__1; ++lls) {
	ls = m - lls + lp1;
/*           ...EXIT */
	if (ls == l) {
	    goto L440;
	}
	test = 0.;
	if (ls != m) {
	    test += (d__1 = e[ls], abs(d__1));
	}
	if (ls != l + 1) {
	    test += (d__1 = e[ls - 1], abs(d__1));
	}
	ztest = test + (d__1 = s[ls], abs(d__1));
	if (ztest != test) {
	    goto L420;
	}
	s[ls] = 0.;
/*           ......EXIT */
	goto L440;
L420:
/* L430: */
	;
    }
L440:
    if (ls != l) {
	goto L450;
    }
    kase = 3;
    goto L470;
L450:
    if (ls != m) {
	goto L460;
    }
    kase = 1;
    goto L470;
L460:
    kase = 2;
    l = ls;
L470:
L480:
    ++l;

/*        PERFORM THE TASK INDICATED BY KASE. */

    switch ((int)kase) {
	case 1:  goto L490;
	case 2:  goto L520;
	case 3:  goto L540;
	case 4:  goto L570;
    }

/*        DEFLATE NEGLIGIBLE S(M). */

L490:
    mm1 = m - 1;
    f = e[m - 1];
    e[m - 1] = 0.;
    i__1 = mm1;
    for (kk = l; kk <= i__1; ++kk) {
	k = mm1 - kk + l;
	t1 = s[k];
	drotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	if (k == l) {
	    goto L500;
	}
	f = -sn * e[k - 1];
	e[k - 1] = cs * e[k - 1];
L500:
	if (wantv) {
	    drot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &
		    cs, &sn);
	}
/* L510: */
    }
    goto L610;

/*        SPLIT AT NEGLIGIBLE S(L). */

L520:
    f = e[l - 1];
    e[l - 1] = 0.;
    i__1 = m;
    for (k = l; k <= i__1; ++k) {
	t1 = s[k];
	drotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	f = -sn * e[k];
	e[k] = cs * e[k];
	if (wantu) {
	    drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L530: */
    }
    goto L610;

/*        PERFORM ONE QR STEP. */

L540:

/*           CALCULATE THE SHIFT. */

/* Computing MAX */
    d__6 = (d__1 = s[m], abs(d__1)), d__7 = (d__2 = s[m - 1], abs(d__2)), 
	    d__6 = max(d__6,d__7), d__7 = (d__3 = e[m - 1], abs(d__3)), d__6 =
	     max(d__6,d__7), d__7 = (d__4 = s[l], abs(d__4)), d__6 = max(d__6,
	    d__7), d__7 = (d__5 = e[l], abs(d__5));
    scale = max(d__6,d__7);
    sm = s[m] / scale;
    smm1 = s[m - 1] / scale;
    emm1 = e[m - 1] / scale;
    sl = s[l] / scale;
    el = e[l] / scale;
/* Computing 2nd power */
    d__1 = emm1;
    b = ((smm1 + sm) * (smm1 - sm) + d__1 * d__1) / 2.;
/* Computing 2nd power */
    d__1 = sm * emm1;
    c = d__1 * d__1;
    shift = 0.;
    if (b == 0. && c == 0.) {
	goto L550;
    }
/* Computing 2nd power */
    d__1 = b;
    shift = sqrt(d__1 * d__1 + c);
    if (b < 0.) {
	shift = -shift;
    }
    shift = c / (b + shift);
L550:
    f = (sl + sm) * (sl - sm) - shift;
    g = sl * el;

/*           CHASE ZEROS. */

    mm1 = m - 1;
    i__1 = mm1;
    for (k = l; k <= i__1; ++k) {
	drotg_(&f, &g, &cs, &sn);
	if (k != l) {
	    e[k - 1] = f;
	}
	f = cs * s[k] + sn * e[k];
	e[k] = cs * e[k] - sn * s[k];
	g = sn * s[k + 1];
	s[k + 1] = cs * s[k + 1];
	if (wantv) {
	    drot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	drotg_(&f, &g, &cs, &sn);
	s[k] = f;
	f = cs * e[k] + sn * s[k + 1];
	s[k + 1] = -sn * e[k] + cs * s[k + 1];
	g = sn * e[k + 1];
	e[k + 1] = cs * e[k + 1];
	if (wantu && k < *n) {
	    drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L560: */
    }
    e[m - 1] = f;
    ++iter;
    goto L610;

/*        CONVERGENCE. */

L570:

/*           MAKE THE SINGULAR VALUE  POSITIVE. */

    if (s[l] >= 0.) {
	goto L580;
    }
    s[l] = -s[l];
    if (wantv) {
	dscal_(p, &c_b44, &v[l * v_dim1 + 1], &c__1);
    }
L580:

/*           ORDER THE SINGULAR VALUE. */

L590:
    if (l == mm) {
	goto L600;
    }
/*           ...EXIT */
    if (s[l] >= s[l + 1]) {
	goto L600;
    }
    t = s[l];
    s[l] = s[l + 1];
    s[l + 1] = t;
    if (wantv && l < *p) {
	dswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
    }
    if (wantu && l < *n) {
	dswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
    }
    ++l;
    goto L590;
L600:
    iter = 0;
    --m;
L610:
    goto L360;
L620:
    return 0;
} /* dsvdc_ */

int drot_(long *n, double *dx, long *incx, 
	double *dy, long *incy, double *c, double *s)
{


    /* System generated locals */
    long i__1;

    /* Local variables */
    static long i;
    static double dtemp;
    static long ix, iy;


/*     applies a plane rotation.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
           to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp = *c * DX(ix) + *s * DY(iy);
	DY(iy) = *c * DY(iy) - *s * DX(ix);
	DX(ix) = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp = *c * DX(i) + *s * DY(i);
	DY(i) = *c * DY(i) - *s * DX(i);
	DX(i) = dtemp;
/* L30: */
    }
    return 0;
} /* drot_ */

int dswap_(long *n, double *dx, long *incx, 
	double *dy, long *incy)
{


    /* System generated locals */
    long i__1;

    /* Local variables */
    static long i, m;
    static double dtemp;
    static long ix, iy, mp1;


/*     interchanges two vectors.   
       uses unrolled loops for increments equal one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
           to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp = DX(ix);
	DX(ix) = DY(iy);
	DY(iy) = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1   


         clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	dtemp = DX(i);
	DX(i) = DY(i);
	DY(i) = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 3) {
	dtemp = DX(i);
	DX(i) = DY(i);
	DY(i) = dtemp;
	dtemp = DX(i + 1);
	DX(i + 1) = DY(i + 1);
	DY(i + 1) = dtemp;
	dtemp = DX(i + 2);
	DX(i + 2) = DY(i + 2);
	DY(i + 2) = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */


int dscal_(long *n, double *da, double *dx, 
	long *incx)
{


    /* System generated locals */
    long i__1, i__2;

    /* Local variables */
    static long i, m, nincx, mp1;


/*     scales a vector by a constant.   
       uses unrolled loops for increment equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	DX(i) = *da * DX(i);
/* L10: */
    }
    return 0;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
	DX(i) = *da * DX(i);
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 5) {
	DX(i) = *da * DX(i);
	DX(i + 1) = *da * DX(i + 1);
	DX(i + 2) = *da * DX(i + 2);
	DX(i + 3) = *da * DX(i + 3);
	DX(i + 4) = *da * DX(i + 4);
/* L50: */
    }
    return 0;
} /* dscal_ */

double dnrm2_(long *n, double *x, long *incx)
{


    /* System generated locals */
    long i__1, i__2;
    double ret_val, d__1;

    /* Local variables */
    static double norm, scale, absxi;
    static long ix;
    static double ssq;


/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DNRM2 := sqrt( x'*x )   



    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   


    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
	norm = 0.;
    } else if (*n == 1) {
	norm = abs(X(1));
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    if (X(ix) != 0.) {
		absxi = (d__1 = X(ix), abs(d__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */

static double c_b44 = -1;

static double c_b4 = 1.;

/* Subroutine */ int drotg_(double *da, double *db, double *c, 
	double *s)
{


    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    static double r, scale, z, roe;


/*     construct givens plane rotation.   
       jack dongarra, linpack, 3/11/78. */


    roe = *db;
    if (abs(*da) > abs(*db)) {
	roe = *da;
    }
    scale = abs(*da) + abs(*db);
    if (scale != 0.) {
	goto L10;
    }
    *c = 1.;
    *s = 0.;
    r = 0.;
    z = 0.;
    goto L20;
L10:
/* Computing 2nd power */
    d__1 = *da / scale;
/* Computing 2nd power */
    d__2 = *db / scale;
    r = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    r = d_sign(&c_b4, &roe) * r;
    *c = *da / r;
    *s = *db / r;
    z = 1.;
    if (abs(*da) > abs(*db)) {
	z = *s;
    }
    if (abs(*db) >= abs(*da) && *c != 0.) {
	z = 1. / *c;
    }
L20:
    *da = r;
    *db = z;
    return 0;
} /* drotg_ */

static long c__9 = 9;
static long c__3 = 3;
static long c__5 = 5;

/* Subroutine */ int ehg182_(long *i)
{
	static string msg100 = " : wrong version number in lowesd.  Probably typo in caller.";
	static string msg101 = " : d>dMAX in ehg131.  Need to recompile with increased dimensions.";
	static string msg102 = " : liv too small.   (Discovered by lowesd)";
	static string msg103 = " : lv too small.    (Discovered by lowesd)";
	static string msg104 = " : alpha too small.  fewer data values than degrees of freedom.";
	static string msg105 = " : k>d2MAX in ehg136.  Need to recompile with increased dimensions.";
	static string msg106 = " : lwork too small";
	static string msg107 = " : invalid value for kernel";
	static string msg108 = " : invalid value for ideg";
	static string msg109 = " : lowstt only applies when kernel=1.";
	static string msg110 = " : not enough extra workspace for robustness calc";
	static string msg120 = " : zero-width neighborhood. make alpha bigger";
	static string msg121 = " : all data on boundary of neighborhood. make alp";
	static string msg122 = " : extrapolation not allowed with blending";
	static string msg123 = " : ihat=1 (diag L) in l2fit only makes sense if z=x (eval=data).";
	static string msg171 = " : lowesd must be called first.";
	static string msg172 = " : lowesf must not come between lowesb and lowese, lowesr, or lowesl.";
	static string msg173 = " : lowesb must come before lowese, lowesr, or lowesl.";
	static string msg174 = " : lowesb need not be called twice.";
	static string msg180 = " : nv>nvmax in cpvert.";
	static string msg181 = " : nt>20 in eval.";
	static string msg182 = " : svddc failed in l2fit.";
	static string msg183 = " : didnt find edge in vleaf.";
	static string msg184 = " : zero-width cell found in vleaf.";
	static string msg185 = " : trouble descending to leaf in vleaf.";
	static string msg186 = " : insufficient workspace for lowesf.";
	static string msg187 = " : insufficient stack space";
	static string msg188 = " : lv too small for computing explicit L";
	static string msg191 = " : computed trace L was negative; something is wrong!";
	static string msg192 = " : computed delta was negative; something is wrong!";
	static string msg193 = " : workspace in loread appears to be corrupted";
	static string msg194 = " : trouble in l2fit/l2tr";
	static string msg195 = " : only constant, linear, or quadratic local models allowed";
	static string msg196 = " : degree must be at least 1 for vertex influence matrix";
	static string msg999 = " : not yet implemented";

    if (*i == 100)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + msg100 );
    if (*i == 101)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + msg101 );
    if (*i == 102)
	throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + msg102 );
    if (*i == 103)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg103 );
    if (*i == 104)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg104 );
    if (*i == 105)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg105 );
    if (*i == 106)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg106 );
    if (*i == 107)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg107 );
    if (*i == 108)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg108 );
    if (*i == 109)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg109 );
    if (*i == 110)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg110 );
    if (*i == 120)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg120 );
    if (*i == 121)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg121 );
    if (*i == 122)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg122 );
    if (*i == 123)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg123 );
    if (*i == 171)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg171 );
    if (*i == 172)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg172 );
    if (*i == 173)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg173 );
    if (*i == 174)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg174 );
    if (*i == 180)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg180 );
    if (*i == 181)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg181 );
    if (*i == 182)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg182 );
    if (*i == 183)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg183 );
    if (*i == 184)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg184 );
    if (*i == 185)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg185 );
    if (*i == 186)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg186 );
    if (*i == 187)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg187 );
    if (*i == 188)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg188 );
    if (*i == 191)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg191 );
    if (*i == 192)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg192 );
    if (*i == 193)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg193 );
    if (*i == 194)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg194 );
    if (*i == 195)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg195 );
    if (*i == 196)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg196 );
    if (*i == 999)
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + msg999 );
    return 0;
} /* ehg182_ */

/* Subroutine */ int ehg183_(char *s, long *i, long *n, long *inc, 
	long s_len)
{
    /* System generated locals */
    long i_dim1, i_offset, i__1;

    /* Builtin functions */
    long s_wsle(cilist *), do_lio(long *, long *, char *, long), 
	    e_wsle(void);

    /* Local variables */
    static long j;

    /* Fortran I/O blocks */
    static cilist io___37 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    i_dim1 = *inc;
    i_offset = i_dim1 + 1;
    i -= i_offset;

    /* Function Body */
    s_wsle(&io___37);
    do_lio(&c__9, &c__1, s, s_len);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	do_lio(&c__3, &c__1, (char *)&i[j * i_dim1 + 1], (long)sizeof(
		long));
    }
    e_wsle();
    return 0;
} /* ehg183_ */

/* Subroutine */ int ehg184_(char *s, double *x, long *n, long *inc,
	 long s_len)
{
    /* System generated locals */
    long x_dim1, x_offset, i__1;

    /* Builtin functions */
    long s_wsle(cilist *), do_lio(long *, long *, char *, long), 
	    e_wsle(void);

    /* Local variables */
    static long j;

    /* Fortran I/O blocks */
    static cilist io___39 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    x_dim1 = *inc;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    s_wsle(&io___39);
    do_lio(&c__9, &c__1, s, s_len);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	do_lio(&c__5, &c__1, (char *)&x[j * x_dim1 + 1], (long)sizeof(
		double));
    }
    e_wsle();
    return 0;
} /* ehg184_ */

/* Subroutine */ int losave_(long *iunit, long *iv, long *liv, 
	long *lv, double *v)
{
    /* Initialized data */

    static long execnt = 0;

    extern /* Subroutine */ int ehg167_(long *, long *, long *, 
	    long *, long *, long *, double *, long *, 
	    double *, double *);

    /* Parameter adjustments */
    --v;
    --iv;

    /* Function Body */
    ++execnt;
    ehg167_(iunit, &iv[2], &iv[4], &iv[5], &iv[6], &iv[14], &v[iv[11]], &iv[
	    iv[7]], &v[iv[12]], &v[iv[13]]);
    return 0;
} /* losave_ */

/* Subroutine */ int ehg167_(long *iunit, long *d, long *vc, long 
	*nc, long *nv, long *nvmax, double *v, long *a, 
	double *xi, double *vval)
{
    /* System generated locals */
    long v_dim1, v_offset, vval_dim1, vval_offset, i__1, i__2;

    /* Builtin functions */
    long s_wsle(cilist *), do_lio(long *, long *, char *, long), 
	    e_wsle(void);

    /* Local variables */
    static long i, j;

    /* Fortran I/O blocks */
    static cilist io___42 = { 0, 0, 0, 0, 0 };
    static cilist io___44 = { 0, 0, 0, 0, 0 };
    static cilist io___46 = { 0, 0, 0, 0, 0 };
    static cilist io___47 = { 0, 0, 0, 0, 0 };
    static cilist io___48 = { 0, 0, 0, 0, 0 };


    /* Parameter adjustments */
    vval_dim1 = *d + 1;
    vval_offset = vval_dim1;
    vval -= vval_offset;
    --xi;
    --a;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;

    /* Function Body */
    io___42.ciunit = *iunit;
    s_wsle(&io___42);
    do_lio(&c__3, &c__1, (char *)&(*d), (long)sizeof(long));
    do_lio(&c__3, &c__1, (char *)&(*nc), (long)sizeof(long));
    do_lio(&c__3, &c__1, (char *)&(*nv), (long)sizeof(long));
    e_wsle();
    i__1 = *d;
    for (i = 1; i <= i__1; ++i) {
/* L10: */
	io___44.ciunit = *iunit;
	s_wsle(&io___44);
	do_lio(&c__5, &c__1, (char *)&v[i * v_dim1 + 1], (long)sizeof(
		double));
	do_lio(&c__5, &c__1, (char *)&v[*vc + i * v_dim1], (long)sizeof(
		double));
	e_wsle();
    }
    j = 0;
    i__1 = *nc;
    for (i = 1; i <= i__1; ++i) {
	if (a[i] != 0) {
	    io___46.ciunit = *iunit;
	    s_wsle(&io___46);
	    do_lio(&c__3, &c__1, (char *)&a[i], (long)sizeof(long));
	    do_lio(&c__5, &c__1, (char *)&xi[i], (long)sizeof(double));
	    e_wsle();
	} else {
	    io___47.ciunit = *iunit;
	    s_wsle(&io___47);
	    do_lio(&c__3, &c__1, (char *)&a[i], (long)sizeof(long));
	    do_lio(&c__3, &c__1, (char *)&j, (long)sizeof(long));
	    e_wsle();
	}
/* L20: */
    }
    i__1 = *nv;
    for (i = 1; i <= i__1; ++i) {
/* L30: */
	io___48.ciunit = *iunit;
	s_wsle(&io___48);
	i__2 = *d;
	for (j = 0; j <= i__2; ++j) {
	    do_lio(&c__5, &c__1, (char *)&vval[j + i * vval_dim1], (long)
		    sizeof(double));
	}
	e_wsle();
    }
    return 0;
} /* ehg167_ */

/* Subroutine */ int lohead_(long *iunit, long *d, long *vc, long 
	*nc, long *nv)
{
    /* Builtin functions */
    long s_rsle(cilist *), do_lio(long *, long *, char *, long), 
	    e_rsle(void), pow_ii(long *, long *);

    /* Fortran I/O blocks */
    static cilist io___49 = { 0, 0, 0, 0, 0 };


    io___49.ciunit = *iunit;
    s_rsle(&io___49);
    do_lio(&c__3, &c__1, (char *)&(*d), (long)sizeof(long));
    do_lio(&c__3, &c__1, (char *)&(*nc), (long)sizeof(long));
    do_lio(&c__3, &c__1, (char *)&(*nv), (long)sizeof(long));
    e_rsle();
    *vc = pow_ii(&c__2, d);
    return 0;
} /* lohead_ */

/* Subroutine */ int loread_(long *iunit, long *d, long *vc, long 
	*nc, long *nv, long *iv, long *liv, long *lv, double *
	v)
{
    /* Initialized data */

    static long execnt = 0;

    extern /* Subroutine */ int ehg182_(long *), ehg168_(long *, 
	    long *, long *, long *, long *, long *, double 
	    *, long *, double *, double *);
    static long bound;

    /* Parameter adjustments */
    --v;
    --iv;

    /* Function Body */
    ++execnt;
    iv[28] = 173;
    iv[2] = *d;
    iv[4] = *vc;
    iv[14] = *nv;
    iv[17] = *nc;
    iv[7] = 50;
    iv[8] = iv[7] + *nc;
    iv[9] = iv[8] + *vc * *nc;
    iv[10] = iv[9] + *nc;
    bound = iv[10] + *nc;
    if (! (bound - 1 <= *liv)) {
	ehg182_(&c__102);
    }
    iv[11] = 50;
    iv[13] = iv[11] + *nv * *d;
    iv[12] = iv[13] + (*d + 1) * *nv;
    bound = iv[12] + *nc;
    if (! (bound - 1 <= *lv)) {
	ehg182_(&c__103);
    }
    ehg168_(iunit, d, vc, nc, nv, nv, &v[iv[11]], &iv[iv[7]], &v[iv[12]], &v[
	    iv[13]]);
    ehg169_(d, vc, nc, nc, nv, nv, &v[iv[11]], &iv[iv[7]], &v[iv[12]], &iv[iv[
	    8]], &iv[iv[9]], &iv[iv[10]]);
    return 0;
} /* loread_ */

/* Subroutine */ int ehg168_(long *iunit, long *d, long *vc, long 
	*nc, long *nv, long *nvmax, double *v, long *a, 
	double *xi, double *vval)
{
    /* System generated locals */
    long v_dim1, v_offset, vval_dim1, vval_offset, i__1, i__2;

    /* Builtin functions */
    long s_rsle(cilist *), do_lio(long *, long *, char *, long), 
	    e_rsle(void);

    /* Local variables */
    static long i, j;

    /* Fortran I/O blocks */
    static cilist io___53 = { 0, 0, 0, 0, 0 };
    static cilist io___54 = { 0, 0, 0, 0, 0 };
    static cilist io___55 = { 0, 0, 0, 0, 0 };


    /* Parameter adjustments */
    vval_dim1 = *d + 1;
    vval_offset = vval_dim1;
    vval -= vval_offset;
    --xi;
    --a;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;

    /* Function Body */
    i__1 = *d;
    for (i = 1; i <= i__1; ++i) {
/* L10: */
	io___53.ciunit = *iunit;
	s_rsle(&io___53);
	do_lio(&c__5, &c__1, (char *)&v[i * v_dim1 + 1], (long)sizeof(
		double));
	do_lio(&c__5, &c__1, (char *)&v[*vc + i * v_dim1], (long)sizeof(
		double));
	e_rsle();
    }
    i__1 = *nc;
    for (i = 1; i <= i__1; ++i) {
/* L20: */
	io___54.ciunit = *iunit;
	s_rsle(&io___54);
	do_lio(&c__3, &c__1, (char *)&a[i], (long)sizeof(long));
	do_lio(&c__5, &c__1, (char *)&xi[i], (long)sizeof(double));
	e_rsle();
    }
    i__1 = *nv;
    for (i = 1; i <= i__1; ++i) {
/* L30: */
	io___55.ciunit = *iunit;
	s_rsle(&io___55);
	i__2 = *d;
	for (j = 0; j <= i__2; ++j) {
	    do_lio(&c__5, &c__1, (char *)&vval[j + i * vval_dim1], (long)
		    sizeof(double));
	}
	e_rsle();
    }
    return 0;
} /* ehg168_ */

/* Subroutine */ int ehg170_(long *k, long *d, long *vc, long *nv,
	 long *nvmax, long *nc, long *ncmax, long *a, long *c, 
	long *hi, long *lo, double *v, double *vval, double 
	*xi)
{
    /* Initialized data */

    static long execnt = 0;

    /* Format strings */
    static char fmt_50[] = "(\002      double precision z(\002,i2,\002)\002)";
    static char fmt_51[] = "(\002      long a(\002,i5,\002), c(\002,i3"
	    ",\002,\002,i5,\002)\002)";
    static char fmt_52[] = "(\002      long hi(\002,i5,\002), lo(\002,i5"
	    ",\002)\002)";
    static char fmt_53[] = "(\002      double precision v(\002,i5,\002,\002,"
	    "i2,\002)\002)";
    static char fmt_54[] = "(\002      double precision vval(0:\002,i2,\002"
	    ",\002,i5,\002)\002)";
    static char fmt_55[] = "(\002      double precision xi(\002,i5,\002)\002)"
	    ;
    static char fmt_56[] = "(\002      double precision ehg128\002)";
    static char fmt_57[] = "(\002      data d,vc,nv,nc /\002,i2,\002,\002,"
	    "i3,\002,\002,i5,\002,\002,i5,\002/\002)";
    static char fmt_58[] = "(\002      data a(\002,i5,\002) /\002,i5,\002"
	    "/\002)";
    static char fmt_59[] = "(\002      data hi(\002,i5,\002),lo(\002,i5,\002"
	    "),xi(\002,i5,\002) /\002,i5,\002,\002,i5,\002,\002,1pe15.6,\002"
	    "/\002)";
    static char fmt_60[] = "(\002      data c(\002,i3,\002,\002,i5,\002) "
	    "/\002,i5,\002/\002)";
    static char fmt_61[] = "(\002      data vval(0,\002,i5,\002) /\002,1pe15"
	    ".6,\002/\002)";
    static char fmt_62[] = "(\002      data v(\002,i5,\002,\002,i2,\002) "
	    "/\002,1pe15.6,\002/\002)";
    static char fmt_63[] = "(\002      data vval(\002,i2,\002,\002,i5,\002"
	    ") /\002,1pe15.6,\002/\002)";

    /* System generated locals */
    long c_dim1, c_offset, v_dim1, v_offset, vval_dim1, vval_offset, i__1, 
	    i__2;

    /* Builtin functions */
    long s_wsle(cilist *), do_lio(long *, long *, char *, long), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(long *, char *, long),
	     e_wsfe(void);

    /* Local variables */
    static long i, j;

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 0, 0, 0, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___60 = { 0, 0, 0, 0, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_51, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_52, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_53, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_54, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_55, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_56, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_57, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_58, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_59, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_61, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_62, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_63, 0 };
    static cilist io___76 = { 0, 0, 0, 0, 0 };
    static cilist io___77 = { 0, 0, 0, 0, 0 };


    /* Parameter adjustments */
    --xi;
    vval_dim1 = *d + 1;
    vval_offset = vval_dim1;
    vval -= vval_offset;
    v_dim1 = *nvmax;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --lo;
    --hi;
    c_dim1 = *vc;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    --a;

    /* Function Body */
    ++execnt;
    io___58.ciunit = *k;
    s_wsle(&io___58);
    do_lio(&c__9, &c__1, "      double precision function loeval(z)", 41L);
    e_wsle();
    io___59.ciunit = *k;
    s_wsfe(&io___59);
    do_fio(&c__1, (char *)&(*d), (long)sizeof(long));
    e_wsfe();
    io___60.ciunit = *k;
    s_wsle(&io___60);
    do_lio(&c__9, &c__1, "      long d,vc,nv,nc", 24L);
    e_wsle();
    io___61.ciunit = *k;
    s_wsfe(&io___61);
    do_fio(&c__1, (char *)&(*nc), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*vc), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*nc), (long)sizeof(long));
    e_wsfe();
    io___62.ciunit = *k;
    s_wsfe(&io___62);
    do_fio(&c__1, (char *)&(*nc), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*nc), (long)sizeof(long));
    e_wsfe();
    io___63.ciunit = *k;
    s_wsfe(&io___63);
    do_fio(&c__1, (char *)&(*nv), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*d), (long)sizeof(long));
    e_wsfe();
    io___64.ciunit = *k;
    s_wsfe(&io___64);
    do_fio(&c__1, (char *)&(*d), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*nv), (long)sizeof(long));
    e_wsfe();
    io___65.ciunit = *k;
    s_wsfe(&io___65);
    do_fio(&c__1, (char *)&(*nc), (long)sizeof(long));
    e_wsfe();
    io___66.ciunit = *k;
    s_wsfe(&io___66);
    e_wsfe();
    io___67.ciunit = *k;
    s_wsfe(&io___67);
    do_fio(&c__1, (char *)&(*d), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*vc), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*nv), (long)sizeof(long));
    do_fio(&c__1, (char *)&(*nc), (long)sizeof(long));
    e_wsfe();
    i__1 = *nc;
    for (i = 1; i <= i__1; ++i) {
	io___69.ciunit = *k;
	s_wsfe(&io___69);
	do_fio(&c__1, (char *)&i, (long)sizeof(long));
	do_fio(&c__1, (char *)&a[i], (long)sizeof(long));
	e_wsfe();
	if (a[i] != 0) {
	    io___70.ciunit = *k;
	    s_wsfe(&io___70);
	    do_fio(&c__1, (char *)&i, (long)sizeof(long));
	    do_fio(&c__1, (char *)&i, (long)sizeof(long));
	    do_fio(&c__1, (char *)&i, (long)sizeof(long));
	    do_fio(&c__1, (char *)&hi[i], (long)sizeof(long));
	    do_fio(&c__1, (char *)&lo[i], (long)sizeof(long));
	    do_fio(&c__1, (char *)&xi[i], (long)sizeof(double));
	    e_wsfe();
	}
	i__2 = *vc;
	for (j = 1; j <= i__2; ++j) {
	    io___72.ciunit = *k;
	    s_wsfe(&io___72);
	    do_fio(&c__1, (char *)&j, (long)sizeof(long));
	    do_fio(&c__1, (char *)&i, (long)sizeof(long));
	    do_fio(&c__1, (char *)&c[j + i * c_dim1], (long)sizeof(long))
		    ;
	    e_wsfe();
/* L4: */
	}
/* L3: */
    }
    i__1 = *nv;
    for (i = 1; i <= i__1; ++i) {
	io___73.ciunit = *k;
	s_wsfe(&io___73);
	do_fio(&c__1, (char *)&i, (long)sizeof(long));
	do_fio(&c__1, (char *)&vval[i * vval_dim1], (long)sizeof(double)
		);
	e_wsfe();
	i__2 = *d;
	for (j = 1; j <= i__2; ++j) {
	    io___74.ciunit = *k;
	    s_wsfe(&io___74);
	    do_fio(&c__1, (char *)&i, (long)sizeof(long));
	    do_fio(&c__1, (char *)&j, (long)sizeof(long));
	    do_fio(&c__1, (char *)&v[i + j * v_dim1], (long)sizeof(
		    double));
	    e_wsfe();
	    io___75.ciunit = *k;
	    s_wsfe(&io___75);
	    do_fio(&c__1, (char *)&j, (long)sizeof(long));
	    do_fio(&c__1, (char *)&i, (long)sizeof(long));
	    do_fio(&c__1, (char *)&vval[j + i * vval_dim1], (long)sizeof(
		    double));
	    e_wsfe();
/* L6: */
	}
/* L5: */
    }
    io___76.ciunit = *k;
    s_wsle(&io___76);
    do_lio(&c__9, &c__1, "      loeval=ehg128(z,d,nc,vc,a,xi,lo,hi,c,v,nv,vv"
	    "al)", 53L);
    e_wsle();
    io___77.ciunit = *k;
    s_wsle(&io___77);
    do_lio(&c__9, &c__1, "      end", 9L);
    e_wsle();
    return 0;
} /* ehg170_ */

CC_END_NAMESPACE()