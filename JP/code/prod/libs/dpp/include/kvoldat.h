/***************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	kvoldat.h
 * Function:	
 * Author:	Christian Daher - David Liu
 ***************************************************************/
#ifndef	_kvoldat_H
#define	_kvoldat_H
#include "kstdinc.h"
#include "ktypes.h"
#include "kzcurve.h"
#include "kmodpar.h"

#include "crxwrpio.h"

/**
 * Enumeration of the basis spread type.
 */
enum KSpd { SUB_SPREAD, PER_SPREAD, ADD_SPREAD };

//--------------------------------------------------------------

/**
 * Class for a calibration index (i.e. 3mCms, 5yFix, etc.)
 * CMS or final, specified either by a date or an interval.
 */

class	KVolCalibIndex {
public:

	/** Default constructor. */
	KVolCalibIndex();

	/** Default destructor. */
	~KVolCalibIndex();

	/**
	 * Converts a calibration index from a string in a DRW format:
	 * nil calibration<br>
	 * CMS intervals (e.g. 1m, 3m, 10yCms, etc.)<br>
	 * Final maturity (diagonal) intervals (e.g. 10yFix, etc.)<br>
	 * Final maturity date (e.g. 1/2/2002)<br>
	 */
	KVolCalibIndex(const char* idxStr);

	/**
	 * Scans a calibration index from a stream in a DRW format
	 * (same format as KVolCalibIndex(const char*)).
	 */
friend	istream& operator>>(istream& is, KVolCalibIndex& object);

	/** Write to a stream. */
friend	ostream& operator<<(ostream& os, const KVolCalibIndex& object);

	/** Equality operator. */
	bool	operator==(const KVolCalibIndex& object) const;

	/**
	 * Returns TRUE if calibration is nil type.
	 */
	bool	IsNil() const
		{ return (mCalibType == 0); }

	/**
	 * Returns TRUE if final calibration, FALSE otherwise.
	 */
	bool	IsFinal() const
		{ return (mCalibType == 2); }

	/**
	 * Returns the final calibration date (or 0 if N/A, i.e.
	 * nil or CMS calibration).
	 */
	TDate	FinalDate() const
		{ return (mFinalAsDate ? mMatDate : 0L); }

	/**
	 * Returns the calibratio interval (or 0 if N/A).
	 */
	KDateInterval	Tenor() const
		{ return mMatInt; }

private:
				/** 0=nil, 1= cms calibration, 2=final */
	int		mCalibType;
				/** FALSE= final as interval, TRUE= as date */
	int		mFinalAsDate;
				/** Maturity calibration as interval */
	KDateInterval	mMatInt;
				/** Maturity calibration as date */
	TDate		mMatDate;

};




//--------------------------------------------------------------
/**
 * Class that contains a generic volatility diagonal
 * (or curve if all tenors are the same) interpolated from the market
 * environment with a  specific calibration index (i.e. 3mCms, 5yFix, etc.)
 */

class	KVolDiag {
public:
	/** Default constructor. */
	KVolDiag();

	/**
	 * Constructor.
	 */
	KVolDiag(
		const KVector(TDate) &volExpDates,
		const KVector(TDate) &volMatDates,
		const KVector(int) &volFreqs,
		const KVector(double) &volRates);

	/**
	 * Constructor.
	 */
	KVolDiag(
		KVolType	      volType,
		const KVector(TDate)  &volExpDates,
		const KVector(TDate)  &volMatDates,
		const KVector(int)    &volFreqs,
		const KVector(double) &volRates);

	/**
	 * Constructs an interpolated diagonal from a DRW environment
	 * and a volatility calibration index.
	 */
	KVolDiag(
		const KVolCalibIndex &calibIdx,	// (I) calibration index
		TDrWrapperData *drWrapData);	// (I) DRW data 


	/**
	 * Construct basis KVolDiag class consistent with the expiration 
	 * dates of the existing IR KVolDiag
 	 */
void	BasisVolDiag(
		const char*	      volType,
		const KVector(TDate)  &volExpDates,
		const KVector(double) &volRates,
		int    		      volFreqs,
		KVolDiag	      &irVolDiag);

	/**
	 * Read IR Vol from MAW voldiag
	 */
void    IRVolDiag(CRX_INPUT  *crxInput);	// (I) MAW data 
    
	/**
	 * Read IR Vol from MAW voldiag
	 */
void    IRVolDiag(BS_INPUT  *bsInput);	// (I) MAW data 

	/**
	 * Read CR Vol from MAW voldiag
	 */
void    CRVolDiag(CRX_INPUT  *crxInput);	// (I) MAW data 
    
	/**
	 * Read SP Vol from MAW voldiag
	 */
void    SPVolDiag(BS_INPUT  *bsInput);	// (I) MAW data 


	/**
	 * Returns the linearly interpolated volatility at the given date.
	 */
	double	VolInterp(TDate date);


	/** Write to a stream. */
friend	ostream& operator<<(ostream& os, const KVolDiag& volDiag);

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Writes the basis spread volatility curve in yacction format.
	 */
ostream& BasisYacctionWrite( ostream& os, int indent=FALSE);


	/** Check for internal consistency. */
	bool	IsEmpty() const
		{ return (mVolRates.size() == 0); }

	/** Check for internal consistency. */
	bool	IsValid() const;

	/** Returns the length. */
	int	Size() const
		{ return mVolRates.size(); }

	/**
	 * Returns true if the two volatility diagonals have
	 * the same set of dates (but may have different volatility
	 * maturities, frequencies, etc.).
	 */
friend	bool	IsConsistent(const KVolDiag& volDiag1,const KVolDiag& volDiag2);

	/**
	 * Set all volatility points to a same value.
	 */
	KVolDiag&	operator=(double argument);

	/**
	 * Assignment operator
	 */
	KVolDiag& 	operator=(const KVolDiag& volDiag);

//private:
				/** Volatility dates */
	KVector(TDate)	mVolDates;
				/** Volatility maturity dates */
	KVector(double)	mVolMats;
				/** Volatility frequencies */
	KVector(int)    mVolFreqs;
				/** Volatility rates */
	KVector(double) mVolRates;
				/** Percentage or basis point volatility */
	KVolType	mVolType;
};




//--------------------------------------------------------------



/**
 * Calibrates a volatility diagonal using vnfm and returns
 * the array of spot volatilities.
 * The argument "nIRDim" is used to possibly use only a subset
 * of the mean reversion parameters up to dimension nIRDim.
 */

void	
DppCalibIRVolDiag(
	const KMrParam &irMrParam,	// (I) full dimension mr info
	const KSmileParam &irSmilePar,	// (I) ir smile info
	int nIRDim,			// (I) ir dimension
	const KVolDiag &volDiag,	// (I) ir vol info
	const KZCurve &diffCurve,	// (I) diffuse zero curve
	KVector(TDate) &volDates,	// (O) ir spot vol dates
	KVector(KVector(double)) &factVolIR);	// (O) ir factor spot vols


void	
DppCalibSpreadVolDiag(
	const KMrParam	&bsMrParam,	// (I) full dimension mr info
	const KSmileParam &bsSmilePar,	// (I) basis smile info
	int		nBSDim,		// (I) basis dimension
	TDate		baseDate,	// (I) base date of ZC
	KSpd		bsType,		// (I) SUB_SPREAD, ADD_SPREAD or PER_SPREAD
	const KZCurve	&refCurve,	// (I) reference zero curve
	const KZCurve	&bsCurve,	// (I) basis zero curve
	const KDayCc	&refDCC,	// (I) reference rate DCC
	const KDayCc	&bsDCC,		// (I) basis rate DCC
	const KVolDiag	&volDiag,	// (I) basis vol info
	KVector(TDate)	&volDates,	// (O) spd spot vol dates
	KVector(KVector(double)) &factVolBasis);// (O) basis factor spot vols



/**
 * Calculate Credit VNFM approximations, given credit BS vols and interest
 * rate spot vols. If crSmilePar.mNumIter<0, then interpret crVolDiag as
 * literal spot vols, rather than BS vols. factVolCurves[0] are the IR spot
 * volatilities. The credit volatilities are returned as factVolCurves[1]
 */
void DppCalibCRSpreadVolDiag(
    const KMrParam           &crMrParam,    /**< (I) CR tree info            */
    const KMrParam           &irMrParam,    /**< (I) IR tree info            */
    double                   ircrCorr,      /**< (I) IR/CR correlation       */
    const KSmileParam        &crSmilePar,   /**< (I) credit smile info       */
    TDate                    baseDate,      /**< (I) base date of ZC         */
    const KZCurve            &refCurve,     /**< (I) IR zero curve           */
    const KZCurve            &crCurve,      /**< (I) Credit zero curve       */
    const KVolDiag           &crVolDiag,    /**< (I) Input credit BS vols    */
    double                   recovery,      /**< (I) Credit recovery rate    */
    KVector(TDate)           &volDates,     /**< (I/O) spot vol dates        */
    KVector(KVector(double)) &factVolCurves /**< (I/O) spot volatilities     */
    );

#endif


