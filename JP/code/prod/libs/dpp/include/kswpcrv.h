/****************************************************************
 * Module:	PenGuin
 * Submodule:	Swap Rates Data Structure
 * File:	kswpcrv.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_kswpcrv_H
#define	_kswpcrv_H

#include "kzcurve.h"

#define	DEF_SWAPRATE_MAX	30




//--------------------------------------------------------------
/**
 * Class for par curve data structure.
 */

class	KPCurve : private KZCurve {
public:
	/** Enumeration for I/O format. */
	enum KPCurveFmt {
		STD,
		BP,			// In basis points
		LIVERATE,		// Liverate text file format
		BVDRWRAP,
		DELTA,			// Delta position
		ZCDRWRAP};


	/** Constructor. */
		KPCurve();

	/** Copy constructor. */
		KPCurve(const KPCurve& swpcrv);

	/** Constructor form TCurve. */
		KPCurve(const TCurve* zcCurve);

	/** Destructor. */
virtual		~KPCurve();

	/** Copy operator. */
	KPCurve& operator=(const KPCurve& swpcrv);

	/** Copy operator. */
	KPCurve& operator=(const TCurve* zcCurve);

	/** Cast into TCurve*. */
		operator TCurve*() { return mZcCurve; }

	/** Returns type name. */
virtual	const char*	TypeName() { return("PCURVE");};

	/** Read from a stream. */
friend	istream& operator>>(istream& is, KPCurve& crv)
			{ crv.Get(is); return is; }

	/** Write to a stream. */
friend	ostream& operator<<(ostream& os, const KPCurve& crv)
			{ crv.Put(os, FALSE); return os; }

	/** Reads the object from a stream. */
virtual	istream& Get(istream& is);

	/** Writes the object to a stream. */
virtual	ostream& Put(ostream& os, int oneLine=FALSE) const;



	/** Returns TRUE if crv1 and crv2 have same dates. */
friend	int	CheckSameType(const KPCurve& crv1, const KPCurve& crv2);

	/** Empties. */
	void	Empty()	{ destroy();}

	/** Returns TRUE if curve is empty. */
	int	IsEmpty() const { return(mZcCurve == NULL);}

	/** Interpolates the value in years. */
	double	InterpYears(double yrs) const;


	/**
	 * Generates a zero curve from a par curve 
	 * using the 5+I method.
	 */
friend	TCurve*	ZCGenNPI(KPCurve& parCrv);

	/**
	 * Generates a zero curve from a par curve
	 * using the 5+I method.
	 */
friend	void	ZCGenNPI(KZCurve& zcCrv, KPCurve& parCrv);

	/** Sets all values. */
	KPCurve&	operator=(double argument);

	/** Adds a constant. */
	KPCurve&	operator+=(double argument);

	/** Multiplies by a constant. */
	KPCurve&	operator*=(double argument);

	/** Adds another curve. */
	KPCurve&	operator+=(const KPCurve& crv);


	/**
	 * Recaculates all effective swap durations
	 * If the argument upDown is TRUE (FALSE) 
	 * the effective duration is based on a two-sided
	 * (one-sided) tweak.
	 */
	KPCurve&	RecalcDuration(int upDown);

	/**
	 * Recalulates all par yields from a zero curve.
	 */
	void		RecalcYields(TCurve *zcCurve);

	/** Computes and returns the 10Y Equivalent. */
	double		Get10YEq();

	/** Resets base date of curve.
	 * If forward is TRUE, forwards the par curve at the specified date.
	 */
	KPCurve&	ShiftDate(TDate baseDate, int forward);

	/** Shifts all yield by an amount value. */
	void		ShiftRates(double value);


	/** Returns the base date.*/
	TDate&		BaseDate() const
				{ return mZcCurve->fBaseDate; }

	/** Returns the number of items. */
	int&		NumItems() const
				{ return mZcCurve->fNumItems; }

	/** Returns the i-th date. */
	TDate&		Date(int i) const
				{ return mZcCurve->fArray[i].fDate; }

	/** Returns the i-th rate. */
	double&		Rate(int i) const
				{ return mZcCurve->fArray[i].fRate; }


	/**
	 * Returns the i-th effective duration
	 * (the duration must have been recalculated using RecalcDuration).
	 */
	double		EffDur(const KDateInterval& intvl, int upDown);


	/**
	 * Returns the i-th effective duration
	 * (the duration must have been recalculated using RecalcDuration).
	 */
	double&		EffDur(int i) const
				{ return mEffDur[i]; }

	/** Day count convention of curve. */
	TDayCount&	DayCount() const
				{ return mZcCurve->fDayCountConv; }

	/** MM denominator. */
	int		MmDen()	 const
				{ return mMmDen; }

	/** Swap frequency. */
	int		Freq() const
				{ return (int) mZcCurve->fBasis; }

	/** Returns the i-th maturity as interval. */
	TDateInterval	Intvl(int i) const;

	/**
	 * Returns the curve as an STL map<KDateInterval,double>.
	 */
	KMap(KDateInterval,double) IntvlMap() const;


	/** Sets the format for I/O. */
	KPCurve&	SetFormat(KPCurveFmt fmt)
				{ mFormat = fmt; return(*this);}


friend	class		KResultData;


private:
					/** private memory management */
	void	create(
			TDate baseDate,
			int size,
			double basis,
			TDayCount dayCountConv);

					/** private memory management */
	void	destroy();

					/** for MM rates */
	int		mMmDen;
					/** Rate type (M,S) (len fNumItems) */
	char		*mRateType;
					/** Swap duration. */
	double		*mEffDur;


};





void	KPCurveRegress(
	int numFact,			// (I) number of factors
	KPCurve *fact,			// (I) array of factors
	KPCurve &weiMat,		// (I) weights (or NULL) 
	KPCurve &beforeMat,		// (I) start (or NULL) 
	KPCurve &afterMat,		// (I) end 
	double *amplitudes,		// (O) amplitudes [0..numFact-1] 
	double *residual);		// (O) residual 


void	KZCurveExportToWrapper(
	const char *fnam,		// (I) file name to export
	const KPCurve &parCurve,	// (I) par curve
	const KZCurve &zcCurve);	// (I) zero curve


#endif	/* _kswpcrv_H */

