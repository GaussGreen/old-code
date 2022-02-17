/***************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	kmodpar.h
 * Function:	
 * Author:	Christian Daher - David Liu
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kmodpar_H
#define	_kmodpar_H
#include "kstdinc.h"
#include "ktypes.h"
#include "kzcurve.h"
#include "crxwrpio.h"


extern	"C" {
#include "dritkwrp.h"           // TDrWrapper
};

//--------------------------------------------------------------
/**
 * Class for multi factor model parameters.
 * It contains the number of factors, mean reversions, weights
 * and correlations.
 */

class	KMrParam {
public:

	/** Default constructor. */
		KMrParam();

	/** Destructor. */
		~KMrParam();

	/** Read from a stream. */
friend	istream& operator>>(istream& is, KMrParam& object);

	/** Write to a stream. */
friend	ostream& operator<<(ostream& os, const KMrParam& object);

	/** Writes in yacction format. */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/** Writes the basis backbone. */
ostream& BasisYacctionWrite( ostream& os, int indent=FALSE);

	/** Writes the credit MAW tree model data. */
ostream& MAWYacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Reads the model parameters from data input stream and
	 * default DRW environment "drWrapData", performing the
	 * default/overwrites as necessary.
	 * The recognized format for the data file is the standard
	 * <PRE>
	 * #Number of factors
	 * <number>
	 * #Betas
	 * <values> or "nil"
	 * #Alphas
	 * <values> or "nil"
	 * #Rhos
	 * <values> or "nil"
	 * #Ppy
	 * <value> or "nil"
	 * </PRE>
	 */
virtual void	ReadDrw(
		istream& is,			// (I) wrapper data stream
		TDrWrapperData *drWrapData);	// (I) wrapper market data 

    /**
     * Reads the IR model parameters from multi-asset wrapper
     * env.
     */
virtual void	IRMAWMrParam(
		istream& is,			// (I) wrapper data stream
		CRX_INPUT *crxInput);	// (I) wrapper market data 

virtual void	IRMAWMrParam(
		istream& is,			// (I) wrapper data stream
		BS_INPUT *bsInput);	    // (I) wrapper market data 

    /**
     * Reads the CR model parameters from multi-asset wrapper
     * env configured for single name credit
     */
virtual void	CRMAWMrParam(
		istream& is,			// (I) wrapper data stream
		CRX_INPUT *crxInput);	// (I) wrapper market data 

virtual void	SPMAWMrParam(
		istream& is,			// (I) wrapper data stream
		BS_INPUT *bsInput);	    // (I) wrapper market data 

	/**
	 * Reads the parameters from a LIL counted array.
	 * The order of the parameters are:
	 *	ppy, number of stddev, smoothing factor, backbone,
	 *	number of factors, betas, weights, and correlations,
	 */
virtual void	ReadLil(const double *mrParamsL);

	/**
	 * Reads the parameters from Magnet interface.
	 * The array of the model parameters are:
	 * number of factors, betas, weights, and correlations, backbone,
	 * The array of the tree parameters include:
	 * ppy, smooth, numStdevCut
	 */
virtual void	ReadMagnet(const KVector(double) &modParams,
			   const KVector(int)    &treeParams);

	/**
	 * Check validity of class members.
	 */
	bool		IsValid();


	/** Sets all values. */
	KMrParam& operator=(double argument);

	/** Multiplies by a constant. */
	KMrParam& operator+=(double argument);

	/** Multiplies by a constant. */
	KMrParam& operator*=(double argument);

	/** Copies. */
	KMrParam& operator=(const KMrParam& m);

	/** Adds. */
	KMrParam& operator+=(const KMrParam& m);

	/**
	 * Returns the alpha-norm of the parameters set
	 */
	double	AlphaNorm() const ;


	/**
	 * Combine two model parameters using a given correlation coefficient.
	 * This routine performs a simple correlation
	 * and can only handle one-dimensional cases (it does not use Vnfm).
	 * The parameters other than mean reversion (mr backbone, ppy, etc.)
	 * are copied from the first argument.
	 */
friend	KMrParam Correlate(
	KMrParam &mrParam1,	// (I)
	KMrParam &mrParam2,	// (I)
	double corr);		// (I)


				/** Tree number of factors. */
	int	mNumFact;
				/** Mean reversions */
	double	mBeta[3];
				/** Factor weights */
	double	mAlpha[3];
				/** Correlations */
	double	mRho[3];

				/** Backbone Q parameter (0=N, 1=LN) */
	double	mBackboneQ;

				/** Number of periods per year */
	int	mPpy;
				/** Smoothing factor for tree */
	int	mSmoothFact;
				/** Mum stdev to cut tree  */
	int	mNumStdevCut;
                /** applies IR/CR corr adjustment to clean spread **/
    int mCorrAdjFlag;
};





//--------------------------------------------------------------
/**
 * Class for skew model parameters.
 * It contains 2-Q coefficents and the number of volatility
 * calibration iterations.
 * WARNING: INTERNALLY, THE "GOOD" CONVENTION HOLDS, I.E.
 * A Q=0 FOR NORMAL, Q=1 FOR LOGNORMAL.
 */

class	KSmileParam {
public:

	/** Default constructor. */
		KSmileParam();

	/** Destructor. */
		~KSmileParam();

	/** Read from a stream. */
friend	istream& operator>>(istream& is, KSmileParam& object);

	/** Write to a stream. */
friend	ostream& operator<<(ostream& os, const KSmileParam& object);

	/** Write to a stream in Yacction format. */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/** Write the basis smile to a stream in Yacction format. */
ostream& BasisYacctionWrite( ostream& os, int indent=FALSE);

	/** Writes the credit MAW tree model data - number CET iterations. */
ostream& MAWYacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Reads the smile parameters in a DRW data file stream is
	 * and performs the overwrites (given a DRW environment
	 * drWrapData). The recognized format for the data file is
	 * L, N or a series of four numbers (q1, q2, forward shift
	 * and Black volatilities iterations). 
	 */
virtual void	ReadDrw(
		istream& is,			// (I) wrapper data stream
		TDrWrapperData *drWrapData);	// (I) wrapper market data 

    /**
     * Reads the IR smile parameters from multi-asset wrapper
     * env 
     */
virtual void	IRMAWSmileParam(
		istream& is,			// (I) wrapper data stream
		CRX_INPUT *crxInput);	// (I) wrapper market data 

virtual void	IRMAWSmileParam(
		istream& is,			// (I) wrapper data stream
		BS_INPUT *bsInput);	    // (I) wrapper market data 
    /**
     * Reads the CR smile parameters from multi-asset wrapper
     * env configured for single name credit
     */
virtual void	CRMAWSmileParam(
		istream& is,			// (I) wrapper data stream
		CRX_INPUT *crxInput);	// (I) wrapper market data 

virtual void	SPMAWSmileParam(
		istream& is,			// (I) wrapper data stream
		BS_INPUT *bsInput);	    // (I) wrapper market data 


	/**
	 * Reads the parameters from a LIL counted array.
	 * The expected order of the parameters is
	 *	q1 (0=LN, 1=N), q2 (0=LN, 1=N),
	 *	forward shift and number of iterations.
	 * WARNING: THE CONVENTION FOR I/O IS 0=LOGNORMAL, 1=NORMAL.
	 */
virtual void	ReadLil(const double *smileParamsL);

	/**
	 * Reads the parameters from Magnet interface.
	 * The array of the smile parameters are:
	 *	q1 (0=LN, 1=N), q2 (0=LN, 1=N),
	 *	forward shift and number of iterations.
	 */
virtual void	ReadMagnet(const KVector(double) &smileParams);



	/**
	 * Check validity of class members.
	 */
	bool		IsValid();


	/** Sets all values. */
	KSmileParam& operator=(double argument);

	/** Multiplies by a constant. */
	KSmileParam& operator+=(double argument);

	/** Multiplies by a constant. */
	KSmileParam& operator*=(double argument);

	/** Copies. */
	KSmileParam& operator=(const KSmileParam& m);

	/** Adds. */
	KSmileParam& operator+=(const KSmileParam& m);

				/** Q1 parameter (0=N, 1=LN). */
	double	mQ1;
				/** Q2 parameter (0=N, 1=LN). */
	double	mQ2;
				/** Forward shift parameter. */
	double	mQF;

	int	mNumIter;
};

const	double	QNORMAL		= 0e0;
const	double	QLOGNORMAL	= 1e0;





#endif


