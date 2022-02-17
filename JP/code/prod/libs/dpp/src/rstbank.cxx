/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "kutilios.h"

#include "krstbank.h"

//---------------------------------------------------------------

KResetBank::KResetBank()
{
	mLastDate = -LONG_MAX;
}



//---------------------------------------------------------------

void
KResetBank::Insert(
	const KRate& rate,	// (I) floating rate
	TDate resetDate,	// (I) reset date
	double value)		// (I) rate value
{

	KPair(KRate, double) rateValue(rate, value);

	// find corresponding date
	mRateResets.insert(KMultimap(TDate, KPair(KRate, double))::value_type(
			resetDate, rateValue));

	mLastDate = MAX(mLastDate, resetDate);
}



//---------------------------------------------------------------

bool	
KResetBank::Get(
	const KRate& rate,	// (I) floating rate
	TDate resetDate,	// (I) reset date
	double *value)		// (O) rate value
{

const	KPair(KRate,double) *rt;

	KRateTable::iterator itD;

	// After last date, don't bother
	if (resetDate > mLastDate) return (FALSE);

	// Check if reset date is there
	//
	switch (mRateResets.count(resetDate))
	{

	    case 0:
		return(false);
		break;

	    case 1:
		itD  = mRateResets.find(resetDate);

		// Check if rate is there
		//
		rt = &((*itD).second);

		if ((*rt).first == rate)
		{
			*value = (*rt).second;

			return(true);
		}
		else
			return(false);
	
		break;

	    default:

		typedef KMultimap(TDate, KPair(KRate,double))::iterator iter;
		KPair(iter, iter) pos;
 
		// pair.first addresses first occurence
		// pair.second addresses position in which 
		// value no longer occurs
		//
		pos = mRateResets.equal_range(resetDate);
		for (; pos.first != pos.second; pos.first++)
		{
			rt = &((*(pos.first)).second);

			if ((*rt).first == rate)
			{
				*value = (*rt).second;

				return(true);
			}
		}	

		return(false);

	}

}




//---------------------------------------------------------------

ostream&
operator<<(ostream& os, const KResetBank& rbank)
{
	TDate		resetDate;
	double		value;
	int		empty = 1;
  const KPair(KRate,double) *rt;
  const KRate		*rate;

	os << "RESET BANK:\n";
	for (KResetBank::KRateTable::const_iterator
			itD = rbank.mRateResets.begin();
	     		itD != rbank.mRateResets.end();
	     		++itD) {
		resetDate = (*itD).first;

		rt = &((*itD).second);
		rate  = &(rt->first);
		value = rt->second;

		os << format("%10s  %12.8f ",
			GtoFormatDate(resetDate),
			value);
		os << *rate;
		os << endl;
		empty = 0;
	}

	if (empty) {
		os << "(empty)" << endl;
	}

	return (os);
}



//---------------------------------------------------------------
//
ostream& 
KResetBank::YacctionWrite( ostream& os, int indent)
{
	TDate		resetDate;
	double		value;
	int		empty = 1;
  const KPair(KRate,double) *rt;
  const KRate		*rate;

	for (KResetBank::KRateTable::const_iterator
			itD = mRateResets.begin();
		itD != mRateResets.end(); ++itD)
	{
		resetDate = (*itD).first;

		rt = &((*itD).second);
		rate  = &(rt->first);
		value = rt->second;

		os << "RESET(";
		os << format("%10s, %10s, %12.8f",
			rate->GetName(),
			GtoFormatDate(resetDate),
			value);
		os << ");" << endl;

	}

	return (os);
}
