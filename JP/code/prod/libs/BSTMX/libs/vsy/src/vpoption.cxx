/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/

#include "vpoption.h"

#include "kutilios.h"		// ios utilities

extern	"C" {
#include "datelist.h"           // GtoNewDateList 
#include "drltime.h"		// DrlTDateIntervalNil() 
#include "drlstr.h"		// DrlStrCmpCaseIndif
};


//-------------------------------------------------------------
// General constructor.
// Strike pay dates are assumed to be the same as
// settlement dates.

KVPOption::KVPOption(
	const char *name,		    // (I) name
	Type type,			    // (I) CALL,PUT,etc
	TBoolean american,		    // (I) Amer/Euro
	const KVector(TDate)& notifDates,   // (I) notif dates
	const KVector(TDate)& settleDates,  // (I) settle dates
	const KVector(double)& strikes,	    // (I) strikes
	const KDateInterval& notDays,	    // (I) # of notifcation days
	const char* discZcName)		    // (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPOption::KVPOption";

	int		idx, numDates;

    try {

	mType = type;
	mAmerican = american;
	mNotifDays = notDays;
	mDiscZcName = String(discZcName);

	// Check length consistency 
	numDates = notifDates.size();
	if (settleDates.size() != numDates)
	    throw KFailure("%s: # settleDates (%d) != # notifDates (%d).\n",
		routine, settleDates.size(), numDates);
	if (strikes.size() != numDates)
	    throw KFailure("%s: # strikes(%d) != # notifDates (%d).\n",
		routine, strikes.size(), numDates);

	for (idx=0; idx<=numDates-1; idx++) {
		mSettleDates.insert(mSettleDates.end(), settleDates[idx]);
		mStrikeDates.insert(mStrikeDates.end(), settleDates[idx]);
		mNotifDates.insert(mNotifDates.end(),   notifDates[idx]);
		mStrikes.insert(mStrikes.end(),         strikes[idx]);

	}
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// General constructor.
// Include the strike pay dates, which can be different
// from settlement dates.

KVPOption::KVPOption(
	const char *name,		    // (I) name
	Type type,			    // (I) CALL,PUT,etc
	TBoolean american,		    // (I) Amer/Euro
	const KVector(TDate)& notifDates,   // (I) notif dates
	const KVector(TDate)& settleDates,  // (I) settle dates
	const KVector(TDate)& strikeDates,  // (I) strike pay dates
	const KVector(double)& strikes,	    // (I) strikes
	const KDateInterval& notDays,       // (I) # of notifcation days
	const char* discZcName)		    // (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPOption::KVPOption";

	int		idx, numDates;

    try {

	mType = type;
	mAmerican = american;
	mNotifDays = notDays;
	mDiscZcName = String(discZcName);

	// Check length consistency 
	numDates = notifDates.size();
	if (settleDates.size() != numDates)
	    throw KFailure("%s: # settleDates (%d) != # notifDates (%d).\n",
		routine, settleDates.size(), numDates);
	if (strikeDates.size() != numDates)
	    throw KFailure("%s: # strikeDates (%d) != # notifDates (%d).\n",
		routine, strikeDates.size(), numDates);
	if (strikes.size() != numDates)
	    throw KFailure("%s: # strikes(%d) != # notifDates (%d).\n",
		routine, strikes.size(), numDates);

	for (idx=0; idx<=numDates-1; idx++) {
		mSettleDates.insert(mSettleDates.end(), settleDates[idx]);
		mStrikeDates.insert(mStrikeDates.end(), strikeDates[idx]);
		mNotifDates.insert(mNotifDates.end(),   notifDates[idx]);
		mStrikes.insert(mStrikes.end(),         strikes[idx]);

	}
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// General constructor.
// Strikes are specified in percentage, where the conversion is
// notionals*(strikePcts - 1)

KVPOption::KVPOption(
	const char *name,		    // (I) name
	Type type,			    // (I) CALL,PUT,etc
	TBoolean american,		    // (I) Amer/Euro
	const KVector(TDate)& notifDates,   // (I) notif dates
	const KVector(TDate)& settleDates,  // (I) settle dates
	const KVector(TDate)& strikeDates,  // (I) strike pay dates
	const KVector(double)& strikePcts,  // (I) strikes
	const KVector(double)& notionals,   // (I) notionals
	const KDateInterval& notDays,	    // (I) # of notifcation days
	const char* discZcName)		    // (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPOption::KVPOption";

	int		idx, numDates;

    try {

	mType = type;
	mAmerican = american;
	mNotifDays = notDays;
	mDiscZcName = String(discZcName);

	// Check length consistency 
	numDates = notifDates.size();
	if (settleDates.size() != numDates)
	    throw KFailure("%s: # settleDates (%d) != # notifDates (%d).\n",
		routine, settleDates.size(), numDates);
	if (strikeDates.size() != numDates)
	    throw KFailure("%s: # strikeDates (%d) != # notifDates (%d).\n",
		routine, strikeDates.size(), numDates);
	if (strikePcts.size() != numDates)
	    throw KFailure("%s: # strikes(%d) != # notifDates (%d).\n",
		routine, strikePcts.size(), numDates);
	if (notionals.size() != numDates)
	    throw KFailure("%s: # notionals(%d) != # notifDates (%d).\n",
		routine, notionals.size(), numDates);

	for (idx=0; idx<=numDates-1; idx++) {
		mSettleDates.insert(mSettleDates.end(), settleDates[idx]);
		mStrikeDates.insert(mStrikeDates.end(), strikeDates[idx]);
		mNotifDates.insert(mNotifDates.end(),   notifDates[idx]);
		mStrikes.insert(mStrikes.end(), 
				notionals[idx]*(strikePcts[idx]-1e0));
	}
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// Convenience constructor to build an exercise table
// from start date to and date with a given frequency.

KVPOption::KVPOption(
	const char *name,	// (I) name
	Type type,		// (I) CALL,PUT,etc
	TBoolean american,	// (I) Amer/Euro
	TDate startDate,	// (I) This date is not included
	TDate matDate,		// (I) maturity date
	TDateInterval freq,	// (I) frequency as interval
	TBoolean stubAtEnd,	// (I) stub
	const KDateInterval&notDays,	// (I) # of notifcation days
	double strike,          // (I) strike
	const char* discZcName) // (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPOption::KVPOption";

	int		idx, numDates;
	TDateList	*dl = NULL;
	TDate		settleDate, notifDate;

    try {

	mType = type;
	mAmerican = american;
	mNotifDays = notDays;
	mDiscZcName = String(discZcName);

	ASSERT_OR_THROW((dl = GtoNewDateList(
		startDate,
		matDate,
		&freq,
		stubAtEnd)) != NULL);
	numDates = dl->fNumItems;

	for (idx=0; idx<=numDates-1; idx++) {
		settleDate = dl->fArray[idx];
		notifDate = settleDate - notDays;

		mSettleDates.insert(mSettleDates.end(), settleDate);
		mStrikeDates.insert(mStrikeDates.end(), settleDate);
		mNotifDates.insert(mNotifDates.end(), notifDate);
		mStrikes.insert(mStrikes.end(), strike);

	}
	GtoFreeDateList(dl);
    }
    catch (KFailure) {
	GtoFreeDateList(dl);
	throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------

istream&
KVPOption::Get(istream& is, int drw)
{

    try {
	if (drw) {
	    int	idx, numDates;

	    mType = KVPOptionType(getString(is, "KVPOption::Get: mType."));

	    mAmerican = getInt(is, "KVPOption::Get: mAmerican.");
	    if (mAmerican)  {
	    	mNotifDays = getKDateInterval(is, "KVPOption::Get: mNotif.");
	    }

	    numDates = getInt(is, "KVPOption::Get: number of dates.");

	    for (idx=0; idx<numDates; idx++) {
		mNotifDates.insert(mNotifDates.end(),
			getTDate(is, "KVPOption::Get: mNotifDates."));

		mSettleDates.insert(mSettleDates.end(),
			getTDate(is, "KVPOption::Get: mSettleDates."));

		mStrikeDates.insert(mStrikeDates.end(),
			getTDate(is, "KVPOption::Get: mStrikeDates."));

		mStrikes.insert(mStrikes.end(),
			getDouble(is, "KVPOption::Get: strikes."));

	    }

	} else {
		throw KFailure("KVPOption::Get: format N/A.\n");
	}


	return(is);
    }
    catch (KFailure) {
	throw KFailure("KVPOption::Get: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPOption::Put(ostream& os, int indent) const
{
	int	idx, numDates;

    try {
	numDates = mNotifDates.size();
	ASSERT_OR_THROW(mSettleDates.size() == numDates);
	ASSERT_OR_THROW(mStrikes.size() == numDates);


	os << "NAME: `" << GetName() << "'" << endl;
	switch (mType) {
	    case CALL:        os << "TYPE: CALL" << endl;        break;
	    case PUT:         os << "TYPE: PUT" << endl;         break;
	    case TIMING_CALL: os << "TYPE: TIMING_CALL" << endl; break;
	    case TIMING_PUT:  os << "TYPE: TIMING_PUT" << endl;  break;
	}

	if (mAmerican)
	{
		os << "Exercise Type: American." << endl;
		os << "NOTIFADJ: " << mNotifDays << endl;
	}
	else
		os << "Exercise Type: European." << endl;


	os << "NUMDATES: " << numDates << endl;

	os << "    NOTIF        SETTLE     STKPAY     STRIKE " << endl;
	for (idx=0; idx<numDates; idx++) {
		os << format(" %10s %10s %10s %18.6f\n",
			DrlTDatePrint(NULL, mNotifDates[idx]),
			DrlTDatePrint(NULL, mSettleDates[idx]),
			DrlTDatePrint(NULL, mStrikeDates[idx]),
			mStrikes[idx]);
	}
	os << "DISCOUNT CURVE: " << mDiscZcName << endl;

	// Print dependencies
	this->KVPAtom::Put(os);

	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPOption::Get: failed.\n");
    }
}



//---------------------------------------------------------------
// Only one underlying is allowed.
//
ostream&
KVPOption::YacctionWrite(ostream& os, int indent)
{
	int	idx, numDates;

    try {

	if (GetWriteFlag())
	{
	    numDates = mNotifDates.size();
	    ASSERT_OR_THROW(mSettleDates.size() == numDates);
	    ASSERT_OR_THROW(mStrikes.size() == numDates);


	    os << GetName() << "=OPTION(" << endl; 
	    switch (mType) {
	        case CALL:        os << "\"CALL\", ";	break;
	        case PUT:         os << "\"PUT\", ";	break;
	        case TIMING_CALL: os << "\"TIMING_CALL\", "; break;
	        case TIMING_PUT:  os << "\"TIMING_PUT\", ";	break;
	    }

	    if (mAmerican)
		os << "TRUE," << endl;
	    else
		os << "FALSE," << endl;

	    os << "{";

	    for (idx=0; idx<numDates; idx++) {
	        os << format(" %10s %10s %10s",
			GtoFormatDate(mNotifDates[idx]),
			GtoFormatDate(mSettleDates[idx]),
			GtoFormatDate(mStrikeDates[idx]));
		os << formatDouble(" %15.3f", mStrikes[idx]) << endl;
	    }

	    if (NumDep() > 1)
		throw KFailure("KVPOption::YacctionWrite: Only "
			"one underlying is allowed in the wrapper.\n");

	    os << "}," << endl
	       << "\t"   << Dep(0)->GetName() << "," << endl
	       << "\t\"" << mDiscZcName << "\");"
	       << endl << endl;

	    WriteDone();
	}
 
	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPOption::YacctionWrite: failed.\n");
    }
}


//---------------------------------------------------------------



KVPOption::Type KVPOptionType(const char *s)
{
static	char	routine[] = "KVPOption::Type";
	KVPOption::Type		type;

    try {
        switch (toupper(s[0])) {
	case 'C':
		type = KVPOption::CALL;
		break;
	case 'P':
		type = KVPOption::PUT;
		break;
	default:
		if (!DrlStrCmpCaseIndif(s, "TIMING_CALL"))
			type = KVPOption::TIMING_CALL;
		else if (!DrlStrCmpCaseIndif(s, "TIMING_PUT"))
			type = KVPOption::TIMING_PUT;
		else 
			throw KFailure("%s: can't scan `%s'.\n", routine, s);
	}

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
    return (type);
}
















