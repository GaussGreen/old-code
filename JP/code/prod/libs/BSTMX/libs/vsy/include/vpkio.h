/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpkio.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vpkio_H
#define	_vpkio_H

#include "vpbase.h"
#include "krate.h"


//--------------------------------------------------------------
/**
 * Class for knock in/out definition. 
 */


class KVPKnockIO : public KVPInstr {
public:

	/**
	 * Convenience constructor given observation freqency.
	 * Single knock rate index.
	 */
	KVPKnockIO(
		const char 	*name,		// (I) name
		const KKnockIO  &ioType,	// (I) Knock in/out type
		const KKnockIO  &ioWindow,	// (I) Pay in/out
		const SharedPointer<KRate> rateIndex,	// (I) Knock in/out rate index
		const KSmooth	&smooth,       	// (I) Node smoothing flag
		TDate           startDate,      // (I) This date is not included
		TDate           matDate,        // (I) maturity date
		const TDateInterval &freq,	// (I) Obs frequency as interval
		TBoolean        stubAtEnd,      // (I) stub
		const KDateInterval   &notDays,	// (I) # of notif days
 
		const KVector(TDate)  &barrierDates,// (I) barrier dates
		const KVector(double) &barrierLo,// (I) Low barries
		const KVector(double) &barrierHi,// (I) High barries
		const KVector(double) &rebates,	// (I) Rebates
		const char* discZcName);	// (I) discount curve

	/**
	 * Convenience constructor for arbitrary observation and 
	 * settlement dates.  Information on effective observation dates
	 * is contained in the rateIndex's spot offset.
	 * Single knock rate index.
	 */
	KVPKnockIO(
		const char *name,		// (I) name
		const KKnockIO &ioType,		// (I) Knock in/out type
		const KKnockIO &ioWindow,	// (I) Pay in/out
		const SharedPointer<KRate> rateIndex, // (I) Knock in/out rate index
		const KSmooth &smooth,       	// (I) Node smoothing flag

		const KVector(TDate)  &obsDates,      // (I) Observ dates
		const KVector(TDate)  &settleDates,   // (I) Settlement dates
		const KVector(double) &barrierLo,     // (I) Low barries
		const KVector(double) &barrierHi,     // (I) High barries
		const KVector(double) &rebates,       // (I) Rebates

		const char*     discZcName);    // (I) discount curve
	/**

	/**
	 * Convenience constructor for arbitrary observation dates
	 * settlement dates.  Effective observation dates are given
	 * explicitly.
	 * Single knock rate index.
	 */
	KVPKnockIO(
		const char *name,		// (I) name
		const KKnockIO &ioType,		// (I) Knock in/out type
		const KKnockIO &ioWindow,	// (I) Pay in/out
		const SharedPointer<KRate> rateIndex, // (I) Knock in/out rate index
		const KSmooth &smooth,       	// (I) Node smoothing flag

		const KVector(TDate)  &obsDates,      // (I) Observ dates
		const KVector(TDate)  &obsEffDates,   // (I) Observ eff dates
		const KVector(TDate)  &settleDates,   // (I) Settlement dates
		const KVector(double) &barrierLo,     // (I) Low barries
		const KVector(double) &barrierHi,     // (I) High barries
		const KVector(double) &rebates,       // (I) Rebates

		const char*     discZcName);    // (I) discount curve

	/**
	 * Convenience constructor for arbitrary knock rate indices and 
	 * observation dates.
	 */
	KVPKnockIO(
		const char *name,		// (I) name
		const KKnockIO &ioType,		// (I) Knock in/out type
		const KKnockIO &ioWindow,	// (I) Pay in/out
		const KSmooth &smooth,       	// (I) Node smoothing flag

		const KVector(TDate)  &obsDates,      // (I) Observ dates
		const KVector(TDate)  &obsEffDates,   // (I) Observ eff dates
		const KVector(TDate)  &settleDates,   // (I) Settlement dates
		const KVector(SharedPointer<KRate>) &rateIndex,// (I) Knock IO rate index
		const KVector(double) &barrierLo,     // (I) Low barries
		const KVector(double) &barrierHi,     // (I) High barries
		const KVector(double) &rebates,       // (I) Rebates

		const char*     discZcName);    // (I) discount curve

	/**
	 * Convenience constructor for NODES or CONTINUOUS observation types..
	 */
	KVPKnockIO(
		const char 	*name,		// (I) name
		const KKnockIO        &ioType,	// (I) Knock in/out type
		const KKnockIO        &ioWindow,	// (I) Pay in/out
		const SharedPointer<KRate> rateIndex,	// (I) Knock in/out rate index
		const KSmooth		&smooth,       	// (I) Node smoothing flag

		const KObsType	obsType,       	// (I) Observation type
		const KDateInterval   &notDays,       // (I) # of notifcation days
 
		const KVector(TDate)  &barrierDates,  // (I) barrier schedule dates
		const KVector(double) &barrierLo,     // (I) Low barries
		const KVector(double) &barrierHi,     // (I) High barries
		const KVector(double) &rebates,       // (I) Rebates

		const char*	discZcName);    // (I) discount curve


	/**
	 * Destructor
	 */
	~KVPKnockIO();


	/**
	 * Include only the events with observation date >= today.
	 */
virtual void	ValidEvents(TDate today);	// today's date

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPKnockIO");}

	/**
	 * Reads the object from a stream.
	 */
virtual	istream& Get(istream& is, int drw=FALSE);

	/**
	 * Writes the object to a stream.
	 */
virtual	ostream& Put(ostream& os, int indent = FALSE) const;

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KVPKnockIO& asset)
		{ asset.Get(is); return (is);}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KVPKnockIO& asset)
		{ asset.Put(os); return (os);}


	/*
	 * Data
	 */
protected:
					/** Knock in, out or none. */
	KKnockIO	mIOType;
					/** Knock in/out window */
	KKnockIO	mIOWindow;

					/** Node smoothing flag */	
	KSmooth		mSmooth;	

					/** Observation type */
	KObsType	mObsType;	
					/** Knock in/out observation dates */
	KVector(TDate)	mObsDates;
					/** Knock in/out settlement dates */
	KVector(TDate)	mSettleDates;
					/** Knock in/out indices */
	KVector(KRateReset*) mRateIndex;
					/** Low barriers */
	KVector(double)	mBarrierLo;
					/** High barriers */
	KVector(double)	mBarrierHi;
					/** Rebates */
	KVector(double)	mRebates;


					/* So can access */
friend	class 	KVPToolKnockIO;
friend  class   KVPToolKnockIO2Idx;

};



#define	BARRIER_TOL	(1e-7)

#endif	/* _vpkio_H */

