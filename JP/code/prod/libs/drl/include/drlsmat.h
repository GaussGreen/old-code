/****************************************************************
 * Module:	DRL
 * Submodule:	SMAT - Swaption Matrix Data Structure
 * File:	drlsmat.h
 * Function:	Swaption matrix routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drlsmat_H
#define	_drlsmat_H
#include "drlstd.h"		/* long, double, etc ... */

#include <stdio.h>

#include "drlts.h"		/* DCurve */

#if defined(DRL_CLIB)

# include "bastypes.h"

typedef	TMatrix2D		DMatrix2D;
typedef	TTable2D		DTable2D;
typedef	TSwaptionMatrix2D	DSwopMat;

#else

/*--------------------------------------------------------------
 * Two-dimensional matrix data structure.
 * <br><br>
 * These types are defined in the header bastypes.h
 * of the C Analytics Library. <br>
 * <i>Important Remark:</i> The structure follows the "matrix" 
 * convention for indices, i.e. line number first. This
 * is opposite to the spreadsheet convention where the column
 * number usually comes first. <br>
 */

typedef struct _DMatrix2D
{
    int      numDim1;        /* Num elems for 1st index */
    int      numDim2;        /* Num elems for 2nd index */
    double **data;           /* Data in matrix */
} DMatrix2D;


/*--------------------------------------------------------------
 * Two-dimensional table data structure.
 * <br><br>
 * These types are defined in the header bastypes.h
 * of the C Analytics Library. <br>
 * <i>Important Remark:</i> The structure follows the "matrix" 
 * convention for indices, i.e. line number first. This
 * is opposite to the spreadsheet convention where the column
 * number usually comes first. <br>
 */

typedef struct _DTable2D
{
    DMatrix2D *matrix;       /* Surface, func of idx1 & idx2 */
    double    *dim1Values;   /* Correspond to 1st idx of matrix */
                             /* Len  = matrix->numDim1  */
    double    *dim2Values;   /* Correspond to 2nd idx of matrix */
                             /* Len  = matrix->numDim2  */
} DTable2D;


/*t-----------------------------------------------------
 * Swaption matrix data structure.
 * <br><br>
 * These types are defined in the header bastypes.h
 * of the C Analytics Library. <br>
 * <i>Important Remark:</i> The structure follows the "matrix" 
 * convention for indices, i.e. line number first. This
 * is opposite to the spreadsheet convention where the column
 * number usually comes first. <br>
 * <i>Time to Expiration:</i> is represented as a year fraction
 * with ACT/365F convention.<br>
 * <i> Time to Maturity:</i> if the matrix is vertical, represents
 * the rate forward maturity as round time interval (i.e. 0.25=3M, etc.).
 * If the matrix is diagonal, represents the year fraction from
 * today to maturity with ACT/365F convention.
 */

typedef struct _DSwopMat
{
    DTable2D  *table;        /* Usually contains swaption vols */
                             /* Dim1 is years to option expiration*/
                             /* Dim2=years from exp to maturity if diag=F*/
                             /* Dim2=years from today to maturity if diag=T*/
    DBoolean   diagonal;     /* Default is FALSE (i.e. vertical matrix)*/
    long       swapPayFreq;  /* Swap payment frequency*/
} DSwopMat;
/*e*/

#endif


#define	TSWAPTION_MATRIX2D_NEXP(ptr)	((ptr)->table->matrix->numDim1)
#define	TSWAPTION_MATRIX2D_NMAT(ptr)	((ptr)->table->matrix->numDim2)
#define	TSWAPTION_MATRIX2D_EXP(ptr, idx)	((ptr)->table->dim1Values[idx])
#define	TSWAPTION_MATRIX2D_MAT(ptr, idx)	((ptr)->table->dim2Values[idx])
#define	TSWAPTION_MATRIX2D_VOL(ptr, idxExp, idxMat) \
		((ptr)->table->matrix->data[idxExp][idxMat])


/*
 *
 */


extern	DSwopMat*	DrlDSwopMatNew(
	DBoolean diagonal,	/* (I) TRUE=diagonal, FALSE=vertical */
	long freq,		/* (I) rate frequency (1,2,4,12) */
	int nExp,		/* (I) # of expirations */
	double *tExp,		/* (I) array of expirations (or NULL) */
	int nMat,		/* (I) # of maturities */
	double *tMat);		/* (I) array of maturities (or NULL) */

extern	DSwopMat*	DrlDSwopMatNewCopy(
				DSwopMat *fromMatrix);
extern	int	DrlDSwopMatFree(DSwopMat *that);
extern	int	DrlDSwopMatIsSameType(
				DSwopMat *that,
				DSwopMat *mat2);

/**
 ** Algebraic operations
 **/
extern	int	DrlDSwopMatOperScalar(DSwopMat *that,
				char *operation, double value);
extern	int	DrlDSwopMatOperMatrix(DSwopMat *that,
				char *operation, DSwopMat *swMat);

extern	int	DrlDSwopMatNumVal(DSwopMat *that,
				char *operation, double *value);

/**
 ** I/O
 **/

#define	TSWAPTION_MATRIX_FMT_STD	((int) 0x00)
#define	TSWAPTION_MATRIX_FMT_NUM	((int) 0x10)
#define	TSWAPTION_MATRIX_FMT_LON	((int) 0x01)
#define	TSWAPTION_MATRIX_FMT_TXT	((int) 0x02)
#define	TSWAPTION_MATRIX_FMT_PRN	((int) 0x11)
#define	TSWAPTION_MATRIX_FMT_CUR	((int) 0x12)
#define	TSWAPTION_MATRIX_FMT_CURK	((int) 0x13)
#define	TSWAPTION_MATRIX_FMT_CURM	((int) 0x14)
#define	TSWAPTION_MATRIX_FMT_LIST	((int) 0x15)


extern	int	DrlDSwopMatFileRead(
	DSwopMat **that,	/* (O) */
	char *fnam,			/* (I) file name */
	int format);			/* (I) TSWAPTION_MATRIX_FMT_XXX */
extern	int	DrlDSwopMatFileWrite(
	DSwopMat *that,	/* (I) */
	char *fnam,			/* (I) file name */
	int format);			/* (I) TSWAPTION_MATRIX_FMT_XXX */
extern	int	DrlDSwopMatFpRead(
	DSwopMat **that,	/* (O) */
	FILE *fp,			/* (I) file pointer */
	int format);			/* (I) TSWAPTION_MATRIX_FMT_XXX */
extern	int	DrlDSwopMatFpWrite(
	DSwopMat *that,	/* (I) */
	FILE *fp,			/* (I) file pointer */
	int format,			/* (I) TSWAPTION_MATRIX_FMT_XXX */
	...);
extern	int	DrlDSwopMatWrapRead(
	DSwopMat **that,	/* (O) new swaption data structure */
	long *typeL,			/* (I) 'L' */
	double *tMatL,			/* (I) 'F' */
	double *tExpL,			/* (I) 'F' */
	double *volL);			/* (I) 'F' */
extern	int	DrlDSwopMatWrapWrite(
	DSwopMat *that,	/* (I) swaptiuon matrix */
	long *typeL,			/* (O) 'L' */
	double *tMatL,			/* (O) 'F' */
	double *tExpL,			/* (O) 'F' */
	double *volL);			/* (O) 'F' */

extern	int	DrlDSwopMatFpReadBVFmt(
	DSwopMat **that,	/* (O) new matrix */
	FILE *fp,			/* (I) file pointer */
	DDate baseDate,			/* (I) reference date */
	int freq);			/* (I) frequency */




extern	int	DrlDSwopMatCreateTimeLine(
	DSwopMat *that,	/* (I) swaption matrix */
	DDate todayDate,		/* (I) reference date */
	int *nDates,			/* (O) # dates */
	DDate **dates);			/* (O) array of dates */

/*
 * Swaption Volatility Interpolation
 */
extern  int DrlGetInterpCoeffs(
        DSwopMat *that,/* (I) swaption matrix */
        double tExp,            /* (I) time to exp */
        double tMat,            /* (I) maturity interpolation time */
        int *ie,                /* (O) exp idx [4] */
        int *im,                /* (O) mat idx [4] */
        double *w);             /* (O) weights [4] */

extern int DrlGetMatInterpCoeffs(
        DSwopMat *that,/* (I) swaption matrix */
        int ie,                 /* (I) exp idx */
        double tMat,            /* (I) mat to interp */
        int *imlo,              /* (O) lo mat idx */
        int *imhi,              /* (O) hi mat idx */
        double *wmlo,           /* (O) lo mat weight */
        double *wmhi);          /* (O) hi mat weight */

extern	int	DrlDSwopMatInterpMatrix(
	DSwopMat *that,	/* (O) output swaption matrix */
	DSwopMat *fromMatrix,	/* (I) input swaption matrix */
        DDate baseDate,			/* (I) base date */
	int adjointFlag);		/* (I) TRUE=rebucket, FALSE=interp */


extern	int	DrlDSwopMatInterpExpMat(
	DSwopMat *that,/* (B) swaption matrix */
	double *value,		/* (B) interpolated value */
	double tExp,		/* (I) time to exp (ACT/365F) */
	double tMat,		/* (I) maturity interpolation time */
	int adjointFlag);	/* (I) TRUE=adjoint, FALSE=direct*/


extern	int	DrlDSwopMatInterpDate(
	DSwopMat *that,/* (B) swaption matrix */
	double *value,		/* (B) interpolated value */
	DDate baseDate,		/* (I) base date */
	DDate expDate,		/* (I) expiration date */
	double tMat,		/* (I) maturity interpolation time */
	int adjointFlag);	/* (I) TRUE=adjoint, FALSE=direct*/


extern	int	DrlDSwopMatInterpMatTime(
	DSwopMat *that,/* (I) swaption matrix */
	DDate baseDate,		/* (I) base date */
	double tExp,		/* (I) time to exp (ACT/365F) */
	int finalFlag,		/* (I) TRUE=final maturity, FALSE=cms */
	DDate matDate,		/* (I) final maturity date (used if final) */
	DInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double *interpMat,	/* (O) matrix mat interp time */
	double *tMat,		/* (O) forward years to maturity (30/360) */
	int *freq);		/* (O) rate frequency (0,1,2,4,12) */

extern	int	DrlDSwopMatInterpValue(
	DSwopMat *that,/* (B) swaption matrix */
	DDate baseDate,		/* (I) base date */
	double tExp,		/* (I) time to exp (ACT/365F) */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	DDate matDate,		/* (I) final maturity date (used if final) */
	DInterval matInt,	/* (I) fwd rate maturity (used if cms) */
	double *value,		/* (B) interpolated volatility */
	double *fwdMat,		/* (I) fwd maturity for interp (or NULL) */
	int adjointFlag);	/* (I) TRUE=rebucket, FALSE=interp */

extern	int	DrlDSwopMatDateInterpValue(
	DSwopMat *that,/* (B) swaption matrix */
	DDate baseDate,		/* (I) base date */
	DDate expDate,		/* (I) expiration date */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	DDate matDate,		/* (I) final maturity date (used if final) */
	DInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double *value,		/* (B) interpolated volatility */
	int adjointFlag);	/* (I) TRUE=rebucket, FALSE=interp */

extern	int	DrlDSwopMatInterpMatrix(
	DSwopMat *that,	/* (B) swaption matrix */
	DSwopMat *interpMat,	/* (B) inpterp swaption matrix */
        DDate baseDate,			/* (I) base date */
	int adjointFlag);		/* (I) TRUE=rebucket, FALSE=interp */

extern	int	DrlDSwopMatInterpVolCurve(
	DSwopMat *that,/* (B) swaption matrix */
	DDate baseDate,		/* (I) base date */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	DDate matDate,		/* (I) final maturity date (used if final) */
	DInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double tMatMin,		/* (I) Minimum muturity                 */
	int adjointFlag,	/* (I) TRUE=rebucket, FALSE=interp */
	int *nVols,		/* (O) # of vol points (can be NULL) */
	DDate **volDates,	/* (O) array of vol exp dates (can be NULL) */
	double **volExp,	/* (O) array of vol exp time (can be NULL) */
	double **volMat,	/* (O) array of vol fwd mat (can be NULL) */
	int **volFreq,		/* (O) array of vol freq (can be NULL) */
	double **volRates,	/* (B) array of vol (can be NULL) */
	DCurve **volCurve);	/* (B) interpolated DCurve (can be NULL) */

#ifdef	_SKIP
/*
 * Miscellaneous functions
 */
extern	DSwopMat*	DrlDSwopMatShift(
	DSwopMat *fromMatrix,	/* (I) swaption matrix */
	DDate baseDate,			/* (I) */
	int diagonal);			/* (I) TRUE=final, FALSE=CMS */


extern	int	DrlDSwopMatColumnNumActiveVols(
	DSwopMat *that,/* (I) swaption matrix */
	int idx,		/* (I) maturity idx */
	int *nActive);		/* (O) # active vols */
#endif


#endif	/* _drlsmat_H */

