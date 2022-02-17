 /****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	krstbank.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_krstbank_H
#define	_krstbank_H

#include "kstdinc.h"
#include "krate.h"


//--------------------------------------------------------------
/**
 * A class to store reset value for floating rates.
 */

class KResetBank {
public:
	/**
	 * Default constructor.
	 */
	KResetBank();

	/**
	 * Returns TRUE is bank is empty (no resets strored).
	 */
	bool	IsEmpty() const
	{ return (mRateResets.size() == 0); }


	/**
	 * Inserts a rate reset in the bank.
	 */
	void	Insert(
		const KRate& rate,	// (I) floating rate
		TDate resetDate,	// (I) reset date
		double value);		// (I) rate value


	/**
	 * Gets a rate reset from the bank.
	 * If no reset is available for rate at the given
	 * reset date, returns FALSE.
	 * Otherwise, sets value and returns TRUE.
	 */
	bool	Get(
		const KRate& rate,	// (I) floating rate
		TDate resetDate,	// (I) reset date
		double *value);		// (I) rate value


	/**
	 * Gets a rate reset from the bank.
	 * If no reset is available for rate at the given
	 * reset date, returns FALSE.
	 * Otherwise, sets value and returns TRUE.
	 */
	bool	Get(
		const KRateReset& reset,// (I) floating rate reset def
		double *value)		// (I) rate value
	{ return Get(reset.Rate(), reset.ResetDate(), value); }



	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KResetBank& rbank);

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);


private:

	typedef	KMultimap(TDate, KPair(KRate,double))	KRateTable;

	/**
	 * Date map of rate resets (a vector for each date).
	 */
	KRateTable		mRateResets;

	/**
	 * Largest date available (for performance).
	 */
	TDate	mLastDate;
};








#endif




