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

#include "drllineq.h"		/* prototype consistency */


/*f---------------------------------------------------------------------
 * Matrix allocation.
 *
 * <br><br>
 * Allocates and returns a <i> double **</i> matrix of size
 * <i> nl x nc</i>  (<i> nl</i> lines and <i> nc</i> columns).
 * Returns <i> NULL</i> if error.
 */

double**
DrlMatrixNew(
	int nl,			/* (I) # of lines */
	int nc)			/* (I) # of columns */
{
static	char	routine[] = "DrlMatrixNew";
	int	status = FAILURE;
	double	**p;
	int	i, j;
	if ((p = DrlDoubleMatrAlloc(0, nl-1, 0, nc-1)) == NULL)
		goto done;

	for(i=0; i<=nl-1; i++)
	for(j=0; j<=nc-1; j++) {
		p[i][j] = 0e0;
	}

	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: %d x %d matrix failed\n",
		routine, nl, nc);
	    return(NULL);
	}
	return(p);
}



/*f---------------------------------------------------------------------
 * Matrix freeing.
 *
 * <br><br>
 * Frees a <i> double **</i> matrix <i> p</i> of size
 * <i> nl x nc</i>
 * (<i> nl</i> lines and <i> nc</i> columns) allocated by <i> DrlMatrixNew</i>.
 * Returns 0 iff OK.
 */

int
DrlMatrixFree(
	double **p,		/* (I) to be freed */
	int nl,			/* (I) # of lines */
	int nc)			/* (I) # of columns */
{
#ifdef	_SKIP
	int	i;
	for(i=0; i<=nl-1; i++) FREE((void*) p[i]);
	FREE((void*) p) ;
	return(nc=0);
#endif
	DrlDoubleMatrFree(p, 0, nl-1, 0, nc-1);
	return(SUCCESS);
}

/*f---------------------------------------------------------------------
 * Matrix printing.
 *
 * <br><br>
 * Prints a <i> double **</i> matrix <i> p</i>
 * <i> nl x nc</i>  (<i> nl</i> lines and <i> nc</i> columns)
 * on the file pointer <i> fp</i>.
 * If <i> fp</i> is NULL, printd on the error log.
 * Returns 0 iff OK.
 */

int
DrlMatrixPrint(
	FILE *fp,		/* (I) file pointer (or NULL) */
	double **p,		/* (I) matrix */
	int nl,			/* (I) # of lines */
	int nc)			/* (I) # of columns */
{
	int	i, j;
#undef	PRINT
#define	PRINT(fmt, arg)	{if (fp == (FILE*)NULL) DrlErrMsg(fmt, arg);\
			else fprintf(fp, fmt, arg);}

	PRINT("NL:%d\n", nl);
	PRINT("NC:%d\n", nc);
	for(i=0; i<=nl-1; i++) {
		for(j=0; j<=nc-1; j++) {
	    	    PRINT("\t%e", p[i][j]);
		}
	        PRINT("\n", 0);
	}

	return(SUCCESS);
#undef	PRINT
}

/*f---------------------------------------------------------------------
 * Matrix : create diagonal matrix.
 *
 * <br><br>
 * Create and returns a <i> double **</i> diagonal square 
 * matrix of size <i> nl x nl</i> with diagonal elements <i> w</i>.
 * Returns <i> NULL</i> if error.
 */

double**
DrlMatrixNewDiag(
	int nl,			/* (I) # of lines */
	double *w)		/* (I) diagonal elements [0..nl-1] */
{
	double	**p;
	int	i, j;

	if ((p = DrlMatrixNew(nl, nl)) == NULL)
		return(NULL);

	for(i=0; i<=nl-1; i++)
	for(j=0; j<=nl-1; j++) {
		p[i][j] = (i == j ? w[i] : 0e0);
	}

	return(p);
}


/*f---------------------------------------------------------------------
 * Matrix copy.
 *
 * <br><br>
 * Copies the
 * <i> double **</i> matrix <i> fromMat</i> of size <i> nl x nc</i>
 * to a matrix <i> toMat</i>.
 */

int
DrlMatrixCopy(
	double ***toMat,	/* (O) new copy matrix */
	double **fromMat,	/* (I) input matrix */
	int nl,			/* (I) # lines */
	int nc)			/* (I) # columns */
{
static	char	routine[] = "DrlMatrixCopy";
	int	i, j;

	if ((*toMat = DrlMatrixNew(nl, nc)) == NULL) {
		DrlErrMsg("%s: matrix allocation failed\n", routine);
		return(-4);
	}

	for(i=0; i<=nl-1; i++)
	for(j=0; j<=nc-1; j++) {
		(*toMat)[i][j] = fromMat[i][j];
	}

	return(0);
}

/*f---------------------------------------------------------------------
 * Matrix transpose.
 *                        
 * <br><br>
 * Transposes the
 * <i> double **</i> matrix <i> mat</i> of size <i> nl x nc</i>.
 */

int
DrlMatrixTranspose(
	double ***mat,		/* (I/O) matrix to transpose */
	int nl,			/* (I) # lines */
	int nc)			/* (I) # columns */
{
static	char	routine[] = "DrlMatrixTranspose";
	register int	i, j;
	double	**tMat;

	if ((tMat = DrlMatrixNew(nc, nl)) == NULL) {
		DrlErrMsg("%s: matrix allocation failed\n", routine);
		return(-4);
	}

	for(i=0; i<=nl-1; i++)
	for(j=0; j<=nc-1; j++) {
		tMat[j][i] = (*mat)[i][j];
	}

	DrlMatrixFree(*mat, nl, nc);
	*mat = tMat;

	return(0);
}


/*f---------------------------------------------------------------------
 * Matrix standard unary operations.
 *
 * <br><br>
 * Performs a matrix operation on
 * <i> double **</i> matrix <i> aMat</i> of size <i> aNl</i>x<i> aNc</i>
 * using as argument the
 * <i> double **</i> matrix <i> bMat</i> of size <i> bNl</i>x<i> bNc</i>.
 * <br> Returns SUCCESS/FAILURE.
 */

int
DrlMatrixUnaryDrlMatrixOper(
	double **aMat,		/* (B) 1st matrix */
	int aNl,		/* (I) # of lines 1st matrix */
	int aNc,		/* (I) # of columns 1st matrix */
	double **bMat,		/* (I) 2nd matrix */
	int bNl,		/* (I) # of lines 2nd matrix */
	int bNc,		/* (I) # of columns 2nd matrix */
	char *oper)		/* (I) operation (=, +=, -=, ) */
{
static	char	routine[] = "DrlMatrixUnaryDrlMatrixOper";
register int	i, j;


	if ((aNl != bNl) || (aNc != bNc)) {
		DrlErrMsg("%s: matrix (%d, %d) not same type as ",
			" matrix (%d, %d)\n",
			routine, aNc, aNl, bNc, bNl);
		return(FAILURE);
	}

	switch (oper[0]) {
	case '+':
		for(i=0; i<=aNl-1; i++)
		for(j=0; j<=aNc-1; j++)
			aMat[i][j] += bMat[i][j];
		break;
	case '-':
		for(i=0; i<=aNl-1; i++)
		for(j=0; j<=aNc-1; j++)
			aMat[i][j] -= bMat[i][j];
		break;
	default:
	    DrlErrMsg("%s: operation `%s' not availble.\n",
		routine, oper);
	    return (FAILURE);
	}

	return(0);
}


/*f---------------------------------------------------------------------
 * Matrix product.
 *
 * <br><br>
 * Computes the product of the
 * <i> double **</i> matrix <i> aMat</i> of size <i> aNl</i>x<i> aNc</i>
 * by the
 * <i> double **</i> matrix <i> bMat</i> of size <i> bNl</i>x<i> bNc</i>.
 * Checks that <i> aNc</i> = <i> bNl</i>.
 * On exit, <i> prodMat</i> points to the new product matrix.
 * Returns 0 iff successful.
 */

int
DrlMatrixProduct(
	double ***prodMat,	/* (O) new product matrix */
	double **aMat,		/* (I) 1st matrix */
	int aNl,		/* (I) # of lines 1st matrix */
	int aNc,		/* (I) # of columns 1st matrix */
	double **bMat,		/* (I) 2nd matrix */
	int bNl,		/* (I) # of lines 2nd matrix */
	int bNc)		/* (I) # of columns 2nd matrix */
{
static	char	routine[] = "DrlMatrixProduct";
	int	i, j, k;
	double	x;


	if (aNc != bNl) {
		DrlErrMsg("%s: can't multiply (aNc = %d, bNl = %d)\n",
			routine, aNc, bNl);
		*prodMat = NULL;
		return(1);
	}

	if ((*prodMat = DrlMatrixNew(aNl, bNc)) == NULL) {
		DrlErrMsg("%s: matrix allocation failed\n", routine);
		return(-4);
	}

	for(i=0; i<=aNl-1; i++)
	for(j=0; j<=bNc-1; j++) {
		x = 0e0;
		for(k=0; k<=aNc-1; k++) {
			x += aMat[i][k]*bMat[k][j];
		}
		(*prodMat)[i][j] = x;
	}

	return(0);
}


/*f---------------------------------------------------------------------
 * Matrix norm.
 *
 * <br><br>
 * Computes and returns the matrix $p$-norm of the
 * <i> double **</i> matrix <i> mat</i> of size
 * <i> nl</i>x<i> nc</i>.
 * The argument <i> p</i> defines the norm: possible chices are
 * <br>
 * <br> p=0 for <br>
 * |A|<sub>infty</sub> = sup {|Av|<sub>infty</sub> / |v|<sub>infty</sub>}
 *   = max<sub>i</sub> sum<sub>j</sub> |a<sub>ij</sub>| ,<br>
 * <br> p=1 for <br>
 * |A|<sub>1</sub> = sup {|Av|<sub>1</sub>/ |v|<sub>1</sub>} 
 *   = max<sub>j</sub> sum<sub>i</sub> |a<sub>ij</sub>| <br>
 * <br> p=2 for <br>
 * |A|<sub>2</sub> = sup {|Av|<sub>2</sub>/ |v|<sub>2</sub>} 
 *   = sqrt{sum <sub>ij</sub> a<sub>ij</sub> <sup>2</sup>}. <br>
 * <br>
 */


double
DrlMatrixNorm(
	double **mat,	/* (I) input matrix */
	int nl,		/* (I) # of lines */
	int nc,		/* (I) # of columns */
	int p)		/* (I) norm type (0,1,2) */
{
static	char	routine[] = "DrlMatrixNorm";
	double		v, w;
	register int	i, j;

	switch (p) {
	case 2:
		v = 0.0;
		for(i=0; i<=nl-1; i++)
		for(j=0; j<=nc-1; j++) 
			v += (mat[i][j])*(mat[i][j]);
		return(sqrt(v));
	case 1:
		v = 0.0;
		for(j=0; j<=nc-1; j++) {
			w = 0.0;
			for(i=0; i<=nl-1; i++)
				w += fabs(mat[i][j]);
			v = MAX(v, w);
		}
		return(v);
	case 0:
		v = 0.0;
		for(i=0; i<=nl-1; i++) {
			w = 0.0;
			for(j=0; j<=nc-1; j++) 
				w += fabs(mat[i][j]);
			v = MAX(v, w);
		}
		break;
	default:
		DrlErrMsg("%s: bad norm type \n", routine);
		v = 1e30;
		break;
	}

	return(v);
}

/*f---------------------------------------------------------------------
 * Matrix L2 scalar product.
 *
 * <br><br>
 * For two <i> nl</i>$x$<i> nc</i> matrices <i> aMat</i> and <i> bMat</i>,
 * computes and returns<br>
 * sum <sub>ij</sub> a<sub>ij</sub> b_<sub>ij</sub>. <br>
 */

double
DrlMatrixL2Product(
	double **aMat,	/* (I) input matrix */
	double **bMat,	/* (I) input matrix */
	int nl,		/* (I) # of lines */
	int nc)		/* (I) # of columns */
{
	double		v;
	register int	i, j;

	v = 0e0;
	for(i=0; i<=nl-1; i++)
	for(j=0; j<=nc-1; j++) 
			v += (aMat[i][j])*(bMat[i][j]);
	return(v);
}


/*f---------------------------------------------------------------------
 * Matrix condition number.
 *
 * <br><br>
 * Computes and returns the value <br>
 * <br>
 * Cond(M) = |M|<sub>p</sub> |M<sup>-1</sup>|<sub>p</sub> <br>
 * <br>
 * where |M|<sub>p</sub> is the $p$-norm of the square matrix <i>M</i>.
 * The argument <i> p</i> defines the norm: possible choices are
 * p=0, p=1 and p=2.
 */

int
DrlMatrixCond(
	double **mat,	/* (I) input square matrix */
	int nl,		/* (I) # of lines */
	int p,		/* (I) norm type */
	double *cond)	/* (O) cond */
{
static	char	routine[] = "DrlMatrixCond";
	double	**invMat;

	if (DrlMatrixInvert(&invMat, mat, nl) != 0) {
		DrlErrMsg("%s: matrix inversion failed\n", routine);
		return(FAILURE);
	}

	*cond  = DrlMatrixNorm(mat, nl, nl,  p)
		* DrlMatrixNorm(invMat, nl, nl, p);

	DrlMatrixFree(invMat, nl, nl);

	return(SUCCESS);
}

/*f---------------------------------------------------------------------
 * Matrix information print.
 *
 * <br><br>
 * Prints the following informations on a square matrix <i> p</i>
 * with <i> nl</i> lines on the file pointer <i> fp</i> (<i> fp</i> is NULL,
 * printd on the error log): the matrix itself, its inverse,
 * its conditionant 1 and 2, its eigen values and eigen vectors.
 * Returns 0 iff OK.
 */

int
DrlMatrixSquareAnalyze(
	FILE *fp,		/* (I) file pointer (or NULL) */
	double **p,		/* (I) matrix */
	int nl)			/* (I) # of lines */
{
	int	status = FAILURE;
	int	i;
	double	**invMat = NULL,
		cond,
		*vap = NULL,
		**vep = NULL;
#undef	PRINT
#define	PRINT(fmt, arg)	{if (fp == (FILE*)NULL) DrlErrMsg(fmt, arg);\
			else fprintf(fp, fmt, arg);}

	/* print matrix */
	PRINT("ANALYZING:\n", 0);
	DrlMatrixPrint(fp, p, nl, nl);


	PRINT("INVERSE:\n", 0);
	if (DrlMatrixInvert(&invMat, p, nl) != 0) 
		goto done;
	DrlMatrixPrint(fp, invMat, nl, nl);


	if (DrlMatrixCond(p, nl, 1, &cond) != SUCCESS)
		goto done;
	PRINT("COND 1: %lf\n", cond);

	if (DrlMatrixCond(p, nl, 2, &cond) != SUCCESS)
		goto done;
	PRINT("COND 2: %lf\n", cond);


	if ((vap = DrlDoubleVectAlloc(0, nl-1)) == NULL)
		goto done;
	if ((vep = DrlMatrixNew(nl, nl)) == NULL)
		goto done;

	if (DrlMatrixRealEigenVect(nl, p, vap, vep) != SUCCESS)
		goto done;
	PRINT("VAP:\n", "");
	for (i=0; i<=nl-1; i++)
		PRINT("\t%e", vap[i]);
	PRINT("\n", 0);
	PRINT("VEP:\n", "");
	DrlMatrixPrint(fp, p, nl, nl);


	status = SUCCESS;
done:
	DrlMatrixFree(invMat, nl, nl);
	DrlDoubleVectFree(vap, 0, nl-1);
	DrlMatrixFree(vep, nl, nl);
	return(status);
}
