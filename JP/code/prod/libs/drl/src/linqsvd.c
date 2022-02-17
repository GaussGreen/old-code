/************************************************************************
 * Module:	DRL - LINEQ
 * Function:	Linear Equations and Matrix Algebra
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <float.h>
#include "drlmem.h"
#include "drlio.h"
#include "drlsort.h"

#include "drllineq.h"		/* prototype consistency */

int
drlSvdcmp(double **a, int m, int n, double w[], double **v);

#define	__DEBUG__
#undef	__DEBUG__

#define	USESVD

/*f---------------------------------------------------------------------
 * Matrix singular value decomposition.
 *
 * <br><br>
 * Performs the singular value decomposition of
 * a <i>m x n</i> matrix <i>A</i>, i.e. 
 * <blockquote>
 * <br> A = U diag(w<sub>0</sub>, ..., w<sub>n</sub>) V<sup>T</sup>.<br>
 * </blockquote>
 * Returns 0 iff successful.
 */

int
DrlMatrixSvdDecomp(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **a,		/* (I) input matrix m x n [0..m-1][0..n-1] */
	double ***u, 		/* (O) output [0..m-1][0..n-1] */
	double **w, 		/* (O) output vap [0..n-1] */
	double ***v) 		/* (O) output [0..n-1][0..n-1] */
{
static	char	routine[] = "DrlMatrixSvdDecomp";
	int	status = FAILURE;


	if (DrlMatrixCopy(u, a, m, n) != SUCCESS)
		goto done;
	if ((*v = DrlMatrixNew(n, n)) == NULL)
		goto done;
	if ((*w = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;


	DrlDoubleMatrCToNr(&a, 1, m, 1, n);
	DrlDoubleMatrCToNr(u, 1, m, 1, n);
	(*w) -= 1;
	DrlDoubleMatrCToNr(v, 1, n, 1, n);

	/*svdcmp(*u, m, n, *w, *v);*/
	drlSvdcmp(*u, m, n, *w, *v);

	DrlDoubleMatrNrToC(&a, 1, m, 1, n);
	DrlDoubleMatrNrToC(u, 1, m, 1, n);
	(*w) += 1;
	DrlDoubleMatrNrToC(v, 1, n, 1, n);


	status = SUCCESS;
	/* made it through OK */
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f---------------------------------------------------------------------
 * Matrix psudo-inversion.
 *
 * <br><br>
 * Computes a pseudo inverse of a <i>m</i> vector 
 * <blockquote>
 * b =(b<sub>0</sub>, ..., b<sub>m</sub>) 
 * </blockquote>
 * from the singular value decomposition of
 * a <i>m x n</i> matrix <i>A</i>, i.e. <br>
 * <blockquote>
 * x = V diag(lambda<sub>1</sub>, ..., lambda<sub>n</sub>) U^T  b. <br>
 * </blockquote>
 * Returns 0 iff successful.
 */

int
DrlMatrixSvdPseudoInverse(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **u, 		/* (I) from SVD [0..m-1][0..n-1] */
	double *lambda,		/* (I) psudo inverse EV [0..n-1] */
	double **v, 		/* (I) from SVD [0..n-1][0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	double *x)		/* (O) solution [0..n-1] */
{
register int	i, j, k;
	/* compute X = V W^-1 tU b */
	for (k=0; k<=n-1; k++) {
	    x[k] = 0e0;
	    for (i=0; i<=n-1; i++) 
	    for (j=0; j<=m-1; j++) 
		x[k] += v[k][i] * lambda[i] * u[j][i] * b[j];
	}

	return(SUCCESS);
}

/*f---------------------------------------------------------------------
 * Matrix : eigenvalues regularisation schemes.
 *
 * <br><br>
 * Computes the pseudo inverse of an array 
 * <blockquote>
 * $(w<sub>1</sub>,...,w<sub>n</sub>) of eigenvectors
 * </blockquote>
 * obtained from an SVD. The arguments <i> regType</i>, <i> alpha</i>
 * and <i> param</i> determine the type and amount
 * of regularization applied.<br>
 *
 * If <i> regType</i> is 0, no regulariztion is performed, i.e.
 * <blockquote>
 *     lambda<sub>i</sub>  = 1/w<sub>i</sub> if w<sub>i</sub> not equal to  0, otherwise 0.
 * </blockquote>
 *
 * If <i> regType</i> is 1, performs Tikhonov regulariztion, i.e.
 * <blockquote>
 *    lambda<sub>i</sub>  = {w<sub>i</sub> / (alpha w<sub>max</sub> <sup>2</sup>)
 *                            + w<sub>i</sub>^2}$$
 *    where w<sub>max</sub>=max|w<sub>i</sub>|.
 * </blockquote>
 * If <i> regType</i> is 2, performs Landweber regulariztion, i.e.
 * <blockquote>
 * \lambda<sub>i</sub>  = {1 - (1- (1 -  param * w<sub>i</sub><sup>2</sup>) <sup>1/alpha</sup>
 *                               / w<sub>i</sub>}.
 * </blockquote>
 * The constant is determined by
 * <blockquote>
 * a = param ( 1 / sup w<sub>i</sub><sup>2</sup> )
 * </blockquote>
 * <br>
 * If <i> regType</i> is 3, performs cutoff regulariztion, i.e. 
 * <blockquote>
 * 	lambda<sub>i</sub>  = 1 / w<sub>i</sub>,
 *      if w<sub>i</sub>&gt;= alpha w<sub>max</sub>, otherwise 0.
 * </blockquote>
 * <br>
 * Returns 0 iff successful.
 */

int
DrlMatrixSvdRegularizeEV(
	int n,			/* (I) # of columns */
	double *w, 		/* (I) EV [0..n-1] */
	double *lambda,		/* (O) pseudo inverse EV [0..n-1] */
	int regType,		/* (I) see below */
	double alpha,		/* (I) see below */
	double param)		/* (I) see below */
{
static	char	routine[] = "DrlMatrixSvdRegularizeEV";
	int	status = FAILURE;
register int	i;
	double	a, evmax;

#ifdef	__DEBUG__
	fprintf(stdout, "%s: regularization %d alpha=%8.4g \n", 
		routine, regType, alpha);
        DrlFPrintDoubleVect(stdout, "\t%8.4g", w, n);
#endif

	switch (regType) {
	case 0:
		/* no regularization */
		for (i=0; i<=n-1; i++) {
			if (IS_ALMOST_ZERO(fabs(w[i])))
				lambda[i] = 0e0;
			else
				lambda[i] = 1e0/w[i];
		}
		break;

	case 1:
		/* compute largest ev (absolute value) */
		evmax = 0e0;
		for (i=0; i<=n-1; i++)
			evmax = MAX(evmax, fabs(w[i]));

		alpha *= evmax*evmax;
		/* Tikhonov */
		for (i=0; i<=n-1; i++) {
			lambda[i] = w[i] / (alpha + w[i]*w[i]);
		}
		break;

	case 2:
		/* Landweber */
		/* compute largest ev (absolute value) */
		evmax = 0e0;
		for (i=0; i<=n-1; i++)
			evmax = MAX(evmax, fabs(w[i]));
		a = param / (evmax * evmax);

		for (i=0; i<=n-1; i++) {
			lambda[i] = 
			    (1e0 - pow(1e0 - a * w[i]*w[i], 1e0/alpha))
				/ w[i];
		}
		break;


	case 3:
		/* cutoff ratio to largest EV */
		/* compute largest ev (absolute value) */
		evmax = 0e0;
		for (i=0; i<=n-1; i++) {
			/*evmax = MAX(fabs(evmax), w[i]);*/
			evmax = MAX(evmax, fabs(w[i]));
		}

		/* cutoff ev <= evmax * ratio */
		for (i=0; i<=n-1; i++) {
			if (fabs(w[i])/evmax <= alpha)
				lambda[i] = 0e0;
			else
				lambda[i] = 1e0/w[i];
		}
		break;

	default:
		DrlErrMsg("%s: bad regularization type %d.\n",
			routine, regType);
		goto done;
	}
#ifdef	__DEBUG__
        DrlFPrintDoubleVect(stdout, "\t%8.4g", lambda, n);
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f---------------------------------------------------------------------
 * Linear system solving with eigenvalue regularisation.
 *
 * <br><br>
 * Solves a possibly ill conditioned linear system
 * <blockquote>
 * sum <sub>j=1,...,n</sub> a<sub>ij</sub> x<sub>j</sub> = b<sub>i</sub>, for i=1,...,m,
 * </blockquote>
 * by singular value decomposition and regularization.
 * See <i> DrlMatrixSvdRegularizeEV</i> for a description
 * of the regularization parameters.
 */

int
DrlRealLinSysSvd(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param)		/* (I) see DrlMatrixSvdRegularizeEV */
{
static	char	routine[] = "DrlRealLinSysSvd";
	int	status = FAILURE;

	double	**u = NULL,	/* [0..m-1][0..n-1] */
		*w = NULL,	/* [0..n-1] */
		**v = NULL;	/* [0..n-1][0..n-1] */


	if (DrlMatrixSvdDecomp(m, n, a, &u, &w, &v) != SUCCESS)
		goto done;

	/* compute pseudo inverse of w */
	if (DrlMatrixSvdRegularizeEV(n, w, w, regType, alpha, param) != SUCCESS)
			goto done;

	/* compute X = V W^-1 tU b */
	if (DrlMatrixSvdPseudoInverse(m, n, u, w, v, b, x) != SUCCESS)
		goto done;

	status = SUCCESS;
	/* made it through OK */
done:
	DrlMatrixFree(u, m, n);
	DrlDoubleVectFree(w, 0, n-1);
	DrlMatrixFree(v, n, n);

	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}




/*----------------------------------------------------------------------
 * NRC2 routines.
 */


#undef	DMAX
#undef	DMIN
#undef	FMAX
#undef	FMIN
#undef	LMAX
#undef	LMIN
#undef	IMAX
#undef	IMIN

#define DMAX(a,b)	((a) > (b) ? (a) : (b))
#define DMIN(a,b)	((a) < (b) ? (a) : (b))
#define FMAX(a,b)	((a) > (b) ? (a) : (b))
#define FMIN(a,b)	((a) < (b) ? (a) : (b))
#define LMAX(a,b)	((a) > (b) ? (a) : (b))
#define LMIN(a,b)	((a) < (b) ? (a) : (b))
#define IMAX(a,b)	((a) > (b) ? (a) : (b))
#define IMIN(a,b)	((a) < (b) ? (a) : (b))


#undef	SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double	sqrarg;
#undef	SQR
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


/*
 *
 */

static	double
drlPythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

/*
 *
 */

int
drlSvdcmp(double **a, int m, int n, double w[], double **v)
{
static	char	routine[] = "drlSvdcmp";
	int status = FAILURE;
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1=NULL;

	rv1=DrlDoubleVectAlloc(1,n);
	if (rv1 == NULL) goto done;

	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=drlPythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) {
				DrlErrMsg("%s: no convergence in 30 "
					"iterations.\n", routine);
				goto done;	
			}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=drlPythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=drlPythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=drlPythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

	status = SUCCESS;
done:
	DrlDoubleVectFree(rv1,1,n);
	if (status != SUCCESS)
		DrlErrMsg("%s: failed.\n", routine);
	return(status);
}



