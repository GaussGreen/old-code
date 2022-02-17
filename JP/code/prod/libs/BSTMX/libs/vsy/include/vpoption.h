/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpoption.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vpoption_H
#define	_vpoption_H

#include "vpbase.h"


//--------------------------------------------------------------
/**
 * Class for option definition. Contains all exercise information
 * (dates, strikes, call/put, American/European, etc.)
 * @author Christian Daher
 * @version 2.0
 */


class KVPOption : public KVPInstr {
public:
	/**
	 *
	 */
	enum Type {
		CALL,		// Call 
		PUT,		// Put
		TIMING_CALL,	// Long Timing Option
		TIMING_PUT	// Short Timing Option
	};


	/**
	 * General constructor.
	 * Strike pay dates are assumed to be the same
	 * as settlement dates.
	 */
	KVPOption(
		const char *name,		    // (I) name
		Type type,			    // (I) CALL,PUT,etc
		TBoolean american,		    // (I) Amer/Euro
		const KVector(TDate)& notifDates,   // (I) notif dates
		const KVector(TDate)& settleDates,  // (I) settle dates
		const KVector(double)& strikes,	    // (I) strikes
		const KDateInterval &notDays,	    // (I) notif interval
		const char* discZcName);	    // (I) discount curve

	/**
	 * General constructor.
	 * Strike pay dates are specified explicitly, which can
	 * be different from settlement dates.
	 */
	KVPOption(
		const char *name,		    // (I) name
		Type type,			    // (I) CALL,PUT,etc
		TBoolean american,		    // (I) Amer/Euro
		const KVector(TDate)& notifDates,   // (I) notif dates
		const KVector(TDate)& settleDates,  // (I) settle dates
		const KVector(TDate)& strikeDates,  // (I) strike pay dates
		const KVector(double)& strikes,	    // (I) strikes
		const KDateInterval &notDays,	    // (I) notif interval
		const char* discZcName);	    // (I) discount curve

	/**
	 * General constructor.
	 * Strike are specified in percentage, where par is given by 100%.
	 * The conversion is notional*(strikePercentage - 1). 
	 */
	KVPOption(
		const char *name,		    // (I) name
		Type type,			    // (I) CALL,PUT,etc
		TBoolean american,		    // (I) Amer/Euro
		const KVector(TDate)& notifDates,   // (I) Notif dates
		const KVector(TDate)& settleDates,  // (I) Settle dates
		const KVector(TDate)& strikeDates,  // (I) Strike pay dates
		const KVector(double)& strikePcts,  // (I) Strikes
		const KVector(double)& notionals,   // (I) Notionals
		const KDateInterval &notDays,	    // (I) Notif interval
		const char* discZcName);	    // (I) Discount curve

	/**
	 * Convenience constructor.
	 */
	KVPOption(
		const char *name,		// (I) name
		Type type,			// (I) call/put
		TBoolean american,		// (I) Amer/Euro
		TDate startDate,		// (I) This date is not included
		TDate matDate,			// (I) maturity date
		TDateInterval freq,		// (I) frequency as interval
		TBoolean stubAtEnd,		// (I) stub
		const KDateInterval &notDays,	// (I) notif interval
		double strike,			// (I) strike
		const char* discZcName);	// (I) discount curve		


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPOption");};

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
friend	istream& operator>>(istream& is, KVPOption& asset)
		{ asset.Get(is); return (is);}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KVPOption& asset)
		{ asset.Put(os); return (os);}

	/**
	 * Reads and return a KVPOption::Type from a string.
	 */
friend	KVPOption::Type KVPOptionType(const char *s);


	/*
	 * Data
	 */
private:
					/** Call/Put/Fwd */
	Type		mType;
					/** American option */
	TBoolean	mAmerican;
					/** Used for amer  */
	KDateInterval	mNotifDays;
					/** Notification dates */
	KVector(TDate)	mNotifDates;
					/** settlement dates */
	KVector(TDate)	mSettleDates;
					/** Strike payment dates. */
	KVector(double)	mStrikeDates;
					/** Strikes */
	KVector(double)	mStrikes;

friend	class	KVPToolOption;
};













#endif	/* _vpoption_H */

