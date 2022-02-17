/****************************************************************
 * Module:	PenGuin
 * File:	kzcurve.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_kzcurve_H
#define	_kzcurve_H

#include "ktypes.h"
#include "kstdinc.h"

#define	DEF_SWAPRATE_MAX	30


//--------------------------------------------------------------
/**
 * Class for zero curve data structure.
 * Wraps the ALIB TCurve structure.
 */

class	KZCurve : public Object {
public:
	/** Enumeration for I/O format. */
	enum KZCurveFmt {
		STD,
		LIVERATE,		// Liverate text file format
		BVDRWRAP,
		ZCDRWRAP};


	/** Constructor. */
		KZCurve();

	/** Copy constructor. */
		KZCurve(const KZCurve& zcrv);

	/**
	 * Constructor from array of dates.
	 */
		KZCurve(
		TDate baseDate,				// (I) 
		const KVector(TDate)& dates,		// (I) 
		const KVector(double)& rates,		// (I) 
		double basis,				// (I) 
		KDayCc dayCountConv,			// (I) 
		int interType = GTO_LINEAR_INTERP);	// (I) 


	/** Constructor from TCurve. */
		KZCurve(const TCurve* zcCurve,
			int interType = GTO_LINEAR_INTERP);

	/** Constructor from TCurve. Shift the base date
	 *  to todayDt, and flat extend the rate to value date. 
	 */
		KZCurve(TCurve* zcCurve, 
			TDate todayDt,
			int interType = GTO_LINEAR_INTERP);


	/** Destructor. */
virtual		~KZCurve();

	/** Copy operator. */
	KZCurve& operator=(const KZCurve& zcrv);

	/** Copy operator. */
	KZCurve& operator=(const TCurve* zcCurve);

	/** Cast into TCurve*. */
		operator TCurve*() const { return mZcCurve; }

	/** Returns type name. */
virtual	const char*	TypeName() { return "ZCURVE";};

	/** Return interpolation type */
virtual int	ZeroInterpType() const { return mZeroInterpType;}

	/** Read from a stream. */
friend	istream& operator>>(istream& is, KZCurve& crv)
			{ crv.Get(is); return is; }

	/** Write to a stream. */
friend	ostream& operator<<(ostream& os, const KZCurve& crv)
			{ crv.Put(os, FALSE); return os; }

	/** Reads the object from a stream. */
virtual	istream& Get(istream& is);

	/** Writes the object to a stream. */
virtual	ostream& Put(ostream& os, int oneLine=FALSE) const;

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE)
	{Put(os, indent); return os;}



	/** Returns TRUE if crv1 and crv2 have same dates. */
friend	int	CheckSameType(const KZCurve& crv1, const KZCurve& crv2);

	/** Returns TRUE if curve is empty. */
	int	IsEmpty() { return(mZcCurve == NULL);}

	/** Interpolates the value in years. */
	double	InterpYears(double yrs) const;

	/** Interpolates and returns the discount factor at a given date. */
	double	DiscFact(TDate discDate) const;

	/** Sets all values. */
	KZCurve&	operator=(double argument);

	/** Adds a constant. */
	KZCurve&	operator+=(double argument);

	/** Multiplies by a constant. */
	KZCurve&	operator*=(double argument);

	/** Adds another curve. */
	KZCurve&	operator+=(const KZCurve& crv);


	/**
	 * Resets base date of curve.
	 * If forward is TRUE, forwards the curve at the specified
	 * date (must be in the future).
	 */
	KZCurve&	ShiftDate(TDate baseDate, int forward);

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

	/** Day count convention of curve. */
	KDayCc		DayCc() const
				{ return KDayCc(mZcCurve->fDayCountConv); }

	/** Swap frequency. */
	double		Basis() const
				{ return (double) mZcCurve->fBasis; }

	/** Swap frequency. */
	int		Freq() const
				{ return (int) mZcCurve->fBasis; }

	/** Returns the i-th maturity as interval. */
	KDateInterval	Intvl(int i) const;

	/** Returns the zero shift between value date
	 *  and today.  Default to 1 if they are the same.
	 */
	double		ZeroShift() const
				{ return mZeroShift; }

	/** Sets the format for I/O. */
	KZCurve&	SetFormat(KZCurveFmt fmt)
				{ mFormat = fmt; return(*this);}



protected:

				/** Used internally. */
	void	create(TDate baseDate, int size, double basis,
				/** Used internally. */
			TDayCount dayCountConv);
	void	destroy();
				/** Used internally. */
	int	isempty() const { return (mZcCurve == NULL); }


					// DATA

					/** Format of the curve. */
	int		mFormat;
					/** Wrapped TCurve structure. */
	TCurve		*mZcCurve;
					/** zero shift between value
					 *  date and today's date
					 */
	double		mZeroShift;

					/** ZC interpolation type */
	int		mZeroInterpType;

};





#endif	/* _kzcurve_H */

