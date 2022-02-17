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
#include "drlstd.h"		/* IntL, FloatL, etc ... */

#include <math.h>
#include <stdio.h>

#include "bastypes.h"

#if !defined(CLIB)

/*t-@CDOC(idxn=TSwaptionMatrix2D,catn=structdef)
 * Two-dimensional matrix data structure.
 * <br><br>
 * These types are defined in the header bastypes.h
 * of the C Analytics Library. <br>
 * <i>Important Remark:</i> The structure follows the "matrix" 
 * convention for indices, i.e. line number first. This
 * is opposite to the spreadsheet convention where the column
 * number usually comes first. <br>
 */

typedef struct _TMatrix2D
{
    int      numDim1;        /* Num elems for 1st index */
    int      numDim2;        /* Num elems for 2nd index */
    double **data;           /* Data in matrix */
} TMatrix2D;

/*e*/


/*t-----------------------------------------------------
 * Two-dimensional table data structure.
 * <br><br>
 * These types are defined in the header bastypes.h
 * of the C Analytics Library. <br>
 * <i>Important Remark:</i> The structure follows the "matrix" 
 * convention for indices, i.e. line number first. This
 * is opposite to the spreadsheet convention where the column
 * number usually comes first. <br>
 */

typedef struct _TTable2D
{
    TMatrix2D *matrix;       /* Surface, func of idx1 & idx2 */
    double    *dim1Values;   /* Correspond to 1st idx of matrix */
                             /* Len  = matrix->numDim1  */
    double    *dim2Values;   /* Correspond to 2nd idx of matrix */
                             /* Len  = matrix->numDim2  */
} TTable2D;
/*e*/


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

typedef struct _TSwaptionMatrix2D
{
    TTablge2D  *table;        /* Usually contains swaption vols */
                             /* Dim1 is years to option expiration*/
                             /* Dim2=years from exp to maturity if diag=F*/
                             /* Dim2=years from today to maturity if diag=T*/
    TBoolean   diagonal;     /* Default is FALSE (i.e. vertical matrix)*/
    long       swapPayFreq;  /* Swap payment frequency*/
} TSwaptionMatrix2D;
/*e*/

#else
# include "bastypes.h"
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


extern	DLL_EXPORT(TSwaptionMatrix2D)	*DrlTSwaptionMatrix2DNew(
	TBoolean diagonal,	/* (I) TRUE=diagonal, FALSE=vertical */
	long freq,		/* (I) rate frequency (1,2,4,12) */
	int nExp,		/* (I) # of expirations */
	double *tExp,		/* (I) array of expirations (or NULL) */
	int nMat,		/* (I) # of maturities */
	double *tMat);		/* (I) array of maturities (or NULL) */

extern	DLL_EXPORT(TSwaptionMatrix2D*)	DrlTSwaptionMatrix2DNewCopy(
				TSwaptionMatrix2D *fromMatrix);
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFree(TSwaptionMatrix2D *that);
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DIsSameType(
				TSwaptionMatrix2D *that,
				TSwaptionMatrix2D *mat2);

/**
 ** Algebraic operations
 **/
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DOperScalar(TSwaptionMatrix2D *that,
				char *operation, double value);
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DOperMatrix(TSwaptionMatrix2D *that,
				char *operation, TSwaptionMatrix2D *swMat);

extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DNumVal(TSwaptionMatrix2D *that,
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


extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFileRead(
	TSwaptionMatrix2D **that,	/* (O) */
	char *fnam,			/* (I) file name */
	int format);			/* (I) TSWAPTION_MATRIX_FMT_XXX */
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFileRead_Mod(
	TSwaptionMatrix2D ***that,	/* (O) */
	char *fnam,			/* (I) file name */
	int format,			/* (I) TSWAPTION_MATRIX_FMT_XXX */
    int nbMatrix);
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFileWrite(
	TSwaptionMatrix2D *that,	/* (I) */
	char *fnam,			/* (I) file name */
	int format);			/* (I) TSWAPTION_MATRIX_FMT_XXX */
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFpRead(
	TSwaptionMatrix2D **that,	/* (O) */
	FILE *fp,			/* (I) file pointer */
	int format);			/* (I) TSWAPTION_MATRIX_FMT_XXX */
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFpRead_Mod(
	TSwaptionMatrix2D ***that,	/* (O) */
	FILE *fp,			/* (I) file pointer */
	int format,			/* (I) TSWAPTION_MATRIX_FMT_XXX */
    int nbMatrix);
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFpWrite(
	TSwaptionMatrix2D *that,	/* (I) */
	FILE *fp,			/* (I) file pointer */
	int format,			/* (I) TSWAPTION_MATRIX_FMT_XXX */
	...);
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DWrapRead(
	TSwaptionMatrix2D **that,	/* (O) new swaption data structure */
	IntL *typeL,			/* (I) 'L' */
	FloatL *tMatL,			/* (I) 'F' */
	FloatL *tExpL,			/* (I) 'F' */
	FloatL *volL);			/* (I) 'F' */
extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DWrapWrite(
	TSwaptionMatrix2D *that,	/* (I) swaptiuon matrix */
	IntL *typeL,			/* (O) 'L' */
	FloatL *tMatL,			/* (O) 'F' */
	FloatL *tExpL,			/* (O) 'F' */
	FloatL *volL);			/* (O) 'F' */

extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DFpReadBVFmt(
	TSwaptionMatrix2D **that,	/* (O) new matrix */
	FILE *fp,			/* (I) file pointer */
	TDate baseDate,			/* (I) reference date */
	int freq);			/* (I) frequency */




extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DCreateTimeLine(
	TSwaptionMatrix2D *that,	/* (I) swaption matrix */
	TDate todayDate,		/* (I) reference date */
	int *nDates,			/* (O) # dates */
	TDate **dates);			/* (O) array of dates */

/*
 * Swaption Volatility Interpolation
 */
extern  DLL_EXPORT(int) DrlGetInterpCoeffs(
        TSwaptionMatrix2D *that,/* (I) swaption matrix */
        double tExp,            /* (I) time to exp */
        double tMat,            /* (I) maturity interpolation time */
        int *ie,                /* (O) exp idx [4] */
        int *im,                /* (O) mat idx [4] */
        double *w);             /* (O) weights [4] */

extern DLL_EXPORT(int) DrlGetMatInterpCoeffs(
        TSwaptionMatrix2D *that,/* (I) swaption matrix */
        int ie,                 /* (I) exp idx */
        double tMat,            /* (I) mat to interp */
        int *imlo,              /* (O) lo mat idx */
        int *imhi,              /* (O) hi mat idx */
        double *wmlo,           /* (O) lo mat weight */
        double *wmhi);          /* (O) hi mat weight */

extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DInterpMatrix(
	TSwaptionMatrix2D *that,	/* (O) output swaption matrix */
	TSwaptionMatrix2D *fromMatrix,	/* (I) input swaption matrix */
        TDate baseDate,			/* (I) base date */
	int adjointFlag);		/* (I) TRUE=rebucket, FALSE=interp */


extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DInterpExpMat(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	double *value,		/* (B) interpolated value */
	double tExp,		/* (I) time to exp (ACT/365F) */
	double tMat,		/* (I) maturity interpolation time */
	int adjointFlag);	/* (I) TRUE=adjoint, FALSE=direct*/


extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DInterpDate(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	double *value,		/* (B) interpolated value */
	TDate baseDate,		/* (I) base date */
	TDate expDate,		/* (I) expiration date */
	double tMat,		/* (I) maturity interpolation time */
	int adjointFlag);	/* (I) TRUE=adjoint, FALSE=direct*/


extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DInterpMatTime(
	TSwaptionMatrix2D *that,/* (I) swaption matrix */
	TDate baseDate,		/* (I) base date */
	double tExp,		/* (I) time to exp (ACT/365F) */
	int finalFlag,		/* (I) TRUE=final maturity, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double *interpMat,	/* (O) matrix mat interp time */
	double *tMat,		/* (O) forward years to maturity (30/360) */
	int *freq);		/* (O) rate frequency (0,1,2,4,12) */

extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DInterpValue(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	TDate baseDate,		/* (I) base date */
	double tExp,		/* (I) time to exp (ACT/365F) */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd rate maturity (used if cms) */
	double *value,		/* (B) interpolated volatility */
	double *fwdMat,		/* (I) fwd maturity for interp (or NULL) */
	int adjointFlag);	/* (I) TRUE=rebucket, FALSE=interp */

extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DDateInterpValue(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	TDate baseDate,		/* (I) base date */
	TDate expDate,		/* (I) expiration date */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double *value,		/* (B) interpolated volatility */
	int adjointFlag);	/* (I) TRUE=rebucket, FALSE=interp */

extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DInterpMatrix(
	TSwaptionMatrix2D *that,	/* (B) swaption matrix */
	TSwaptionMatrix2D *interpMat,	/* (B) inpterp swaption matrix */
        TDate baseDate,			/* (I) base date */
	int adjointFlag);		/* (I) TRUE=rebucket, FALSE=interp */

extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DInterpVolCurve(
	TSwaptionMatrix2D *that,/* (B) swaption matrix */
	TDate baseDate,		/* (I) base date */
	int finalFlag,		/* (I) TRUE=final, FALSE=cms */
	TDate matDate,		/* (I) final maturity date (used if final) */
	TDateInterval matInt,	/* (I) fwd maturity interval (used if cms) */
	double tMatMin,		/* (I) Minimum muturity                 */
	int adjointFlag,	/* (I) TRUE=rebucket, FALSE=interp */
	int *nVols,		/* (O) # of vol points (can be NULL) */
	TDate **volDates,	/* (O) array of vol exp dates (can be NULL) */
	double **volExp,	/* (O) array of vol exp time (can be NULL) */
	double **volMat,	/* (O) array of vol fwd mat (can be NULL) */
	int **volFreq,		/* (O) array of vol freq (can be NULL) */
	double **volRates,	/* (B) array of vol (can be NULL) */
	TCurve **volCurve);	/* (B) interpolated TCurve (can be NULL) */

#ifdef	_SKIP
/*
 * Miscellaneous functions
 */
extern	DLL_EXPORT(TSwaptionMatrix2D*)	DrlTSwaptionMatrix2DShift(
	TSwaptionMatrix2D *fromMatrix,	/* (I) swaption matrix */
	TDate baseDate,			/* (I) */
	int diagonal);			/* (I) TRUE=final, FALSE=CMS */


extern	DLL_EXPORT(int)	DrlTSwaptionMatrix2DColumnNumActiveVols(
	TSwaptionMatrix2D *that,/* (I) swaption matrix */
	int idx,		/* (I) maturity idx */
	int *nActive);		/* (O) # active vols */
#endif


#endif	/* _drlsmat_H */

