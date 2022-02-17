/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	kemmatr.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_kemmatr_H
#define	_kemmatr_H
#include "ktypes.h"

/*
extern "C" {
#include "drlsmat.h"
}
*/

//----------------------------------------------------------------------
/**
 * Expiration/maturity (swaption) matrix class.
 * Wraps the ALIB type TSwaptionMatrix2D.
 */


class	KEmMatrix {
public:
	/** Enumeration for I/O format. */
	enum KEmMatrixFmt {
		DEF,			// Standard CD format
		CUR,
		CURK,
		CURM,
		PERCENT,		// Same but in %
		TXT,
		PRN,
		LIST,
		LIVERATE,		// Liverate format
		BVMAT,			// base vol matrix format for Kap feed
		DRW2BVCURVE,		// DRW II base vol curve.
		DRW2BVMAT,		// DRW II base vol matrix.
		DRW2SWMAT		// DRW II swaption vol.
	};


	/** Returns the type name */
virtual	const char*	TypeName() const { return("SWOPMAT");};

	/** Default constructor. */
		KEmMatrix();

	/** Copy constructor. */
		KEmMatrix(const KEmMatrix& argument);

	/** Copy constructor. */
		KEmMatrix(TSwaptionMatrix2D* argument);

	/** Destructor. */
		~KEmMatrix();

	/**
	 * Checks that object is valid and returns TRUE/FALSE
	 * (validity consists in having increasing expirations
	 * and maturities).
	 * Throws an exception on an invalid object.
	 */
	void	CheckValid();

	/**
	 * Check if matrix is empty.
	 */
	int	IsEmpty() 
		{ return (mMatrix == NULL); }

	/**
	 * Copy operator.
	 */
	KEmMatrix& operator=(const KEmMatrix& argument);

	/**
	 * Copy operator.
	 */
	KEmMatrix& operator=(TSwaptionMatrix2D* argument);

	/**
	 * Copy operator.
	 */
	KEmMatrix& operator=(TCurve* curve);



	/** Read from stream. */
friend	istream& operator>>(istream& is, KEmMatrix& mat)
		{ mat.Get(is); return(is);}

	/** Write to stream. */
friend	ostream& operator<<(ostream& os, const KEmMatrix& mat)
		{ mat.Put(os, FALSE); return(os);}

	/** Reads the object from a stream. */
virtual	istream& Get(istream& is);


	/** Writes the object to a stream. */
virtual	ostream& Put(ostream& os, int oneLine=FALSE) const;


	/** Casting to TSwaptionMatrix2D&. */
	operator	TSwaptionMatrix2D&() { return(*mMatrix);}

	/** Casting to TSwaptionMatrix2D*. */
	operator	TSwaptionMatrix2D*() { return(mMatrix);}


	/**
	 * Sets all values.
	 */
	KEmMatrix& operator=(double argument);

	/**
	 * Adds a constant.
	 */
	KEmMatrix& operator+=(double argument);

	/**
	 * Multiplies by a constant.
	 */
	KEmMatrix& operator*=(double argument);

	/**
	 * Adds another matrix.
	 */
	KEmMatrix& operator+=(const KEmMatrix& crv);

	/**
	 * Subtracts another matrix.
	 */
	KEmMatrix& operator-=(const KEmMatrix& crv);


	/**
	 * Computes and returns the sum.
	 */
	double	Sum() const;

	/**
	 * Computes and returns the sumproduct of the matrix
	 * with another one.
	 */
	double	SumProd(const KEmMatrix& mat) const;


	/** Interpolates the matrix */
	KEmMatrix& Interp(KEmMatrix& mat, TDate baseDate, int interpFlag);

	/** Number of expirations. */
	int&	NExp() const	{ return (mMatrix->table->matrix->numDim1);}

	/** Number of maturities. */
	int&	NMat() const	{ return (mMatrix->table->matrix->numDim2);}

	/** Expirations. */
	double&	TExp(int idxExp) const
			    { return (mMatrix->table->dim1Values[idxExp]);}

	/** Maturities. */
	double&	TMat(int idxMat) const
			    { return (mMatrix->table->dim2Values[idxMat]);}

	/** Volatilities. */
	double&	Val(int idxExp, int idxMat)  const { return
			    (mMatrix->table->matrix->data[idxExp][idxMat]);}


	/** Swap volatility frequency. */
	long&	Freq() const { return (mMatrix->swapPayFreq); }

	/** Vertical (0) or diagonal(1). */
	int&	Diag() const { return (mMatrix->diagonal); }


	/**
	 * Returns TRUE if mat1 and mat2 have same expirations
	 * and maturities.
 	 */
friend	int	CheckSameType(const KEmMatrix& mat1, const KEmMatrix& mat2);

	/** Sets the format for I/O */
	void	SetFormat(KEmMatrixFmt format)
			{ mFormat = format; }

	/**
	 * Interpolates at given expiration/maturity in years.
	 */
	double	Interp(double expYrs, double matYrs);


	/**
	 * Interpolates a matrix for constant or final maturity
	 */
	void 	InterpCurve(
		TDate baseDate,		// (I) base date
		int finalFlag,		// (I) TRUE=final, FALSE=cms 
		int finalDate,		// (I) TRUE=use date, FALSE=interval
		TDate matDate,		// (I) final mat (used if final) 
		KDateInterval matInt,	// (I) fwd mat (used if cms) 
		KDateInterval minMat,	// (I) minimum maturity
		KVector(TDate) &vDates,	// (O) array of vol exp dates
		KVector(double) &vExp,	// (O) array of vol exp time
		KVector(double) &vMat,	// (O) array of vol fwd mat
		KVector(int) &vFreq,	// (O) array of vol freq
		KVector(double) &vRates);// (B) array of vol 

	/**
	 * Computes and returns the maturity interpolation
	 * indices and weights.
	 */
	void	InterpMatCoeffs(
		int ie,			// (I) expiration idx
		double tMat,		// (I) maturoty to interpolate
		int *imlo,		// (O) lo mat idx
		int *imhi,		// (O) hi mat idx
		double *wmlo,		// (O) lo mat weight
		double *wmhi);		// (O) hi mat weight

private:
	KEmMatrixFmt		mFormat;
	TSwaptionMatrix2D	*mMatrix;

	int	cleanup();
	int	isempty() const
			{ return (mMatrix == NULL);}

};

//
//
//

void	KEmMatrixRegress(
	int numSwfact,			// (I) number of vol scen 
	KEmMatrix *swfact,		// (I) swaption scenarios 
	KEmMatrix &weiMat,		// (I) weights matr (or NULL) 
	KEmMatrix &beforeMat,		// (I) start matrix (or NULL) 
	KEmMatrix &afterMat,		// (I) end matrix 
	double *amplitudes,		// (O) amplitudes [0..numSwfact-1] 
	double *residual);		// (O) residual 




/*----------------------------------------------------------------------
 * LEGACY: Please do not use TSwaptionMatrix2D* !
 */

/**@#-*/

typedef	TSwaptionMatrix2D	TSwoptMat;


inline	int	NExp(TSwoptMat *that)
			{ return (that->table->matrix->numDim1);}
inline	int	NMat(TSwoptMat *that)
			{ return (that->table->matrix->numDim2);}
inline	int	Freq(TSwoptMat *that)
			{ return (that->swapPayFreq); }
inline	int	Diag(TSwoptMat *that)
			{ return (that->diagonal); }
inline	double&	TExp(TSwoptMat *that, int idxExp)
			{ return (that->table->dim1Values[idxExp]);}
inline	double&	TMat(TSwoptMat *that, int idxMat)
			{ return (that->table->dim2Values[idxMat]);}
inline	double&	Vol(TSwoptMat *that, int idxExp, int idxMat)
			{ return (that->table->matrix->data[idxExp][idxMat]);}

/**@#+*/


#endif	/* _kemmatr_H */

