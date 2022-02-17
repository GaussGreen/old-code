/************************************************************************
 * Module:	PenGuin
 * File:	vpdefprotect.cxx
 * Function:    
 * Author:      David Liu
 ************************************************************************/

#include "vpdefprotect.h"

#include "kutilios.h"		// ios utilities

extern	"C" {
#include "datelist.h"           // GtoNewDateList 
#include "drltime.h"		// DrlTDateIntervalNil() 
#include "drlstr.h"		// DrlStrCmpCaseIndif
};


//-------------------------------------------------------------
// Convenience constructor to build a protection schedule
// from start date to and date with a given frequency.

KVPDefProtect::KVPDefProtect(
        const char          *name,          // (I) name
        TDate               startDate,      // (I) start date of period
        TDate               endDate,        // (I) end date of period
        TDateInterval       &freq,          // (I) Obs frequency as interval
        DefType             type,           // (I) default protection type
        double              recovery,       // (I) Recovery rate
        const char*         discZcName)     // (I) discount curve
	: KVPInstr(name)
{
static	char	routine[] = "KVPDefProtect::KVPDefProtect";

	int		idx, numDates;
	TDateList	*dl = NULL;
	TDate		stDate, settleDate;
        TBoolean        stubAtEnd = FALSE; 

    try {

	mDiscZcName = String(discZcName);
        mDefType    = type;
        mRecovery   = recovery;
	

        // Needs at least one protection period
        //
        if (startDate == endDate)
            throw KFailure("%s: invalid protection period (start date (%s) "
                           "= end date (%s).\n",
                           routine,
                           GtoFormatDate(startDate), 
                           GtoFormatDate(endDate));

	ASSERT_OR_THROW((dl = GtoNewDateList(
		startDate,
		endDate,
		&freq,
		stubAtEnd)) != NULL);
	numDates = dl->fNumItems;

        //
        // Total numDates-1 periods 
        //
	for (idx=0; idx<numDates-1; idx++) {
		stDate  = dl->fArray[idx];
		settleDate = dl->fArray[idx+1];

		mStartDates.insert(mStartDates.end(), stDate);
		mEndDates.insert(mEndDates.end(), settleDate);
		mSettleDates.insert(mSettleDates.end(), settleDate);
		mRebates.insert(mRebates.end(), 0e0);
	}

	GtoFreeDateList(dl);
    }
    catch (KFailure) {
	GtoFreeDateList(dl);
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// General constructor.

KVPDefProtect::KVPDefProtect(
        const char            *name,         // (I) name
        const KVector(TDate)  &startDates,   // (I) Start dates
        const KVector(TDate)  &endDates,     // (I) End dates
        const KVector(TDate)  &settleDates,  // (I) Settlement dates
        DefType               type,          // (I) default protection type
        double                recovery,      // (I) Recovery rate
        const char*           discZcName)    // (I) discount curve
	: KVPInstr(name)
{
static	char	routine[] = "KVPDefProtect::KVPDefProtect";

	int	idx, numDates;

    try {

	mDiscZcName = String(discZcName);
        mDefType    = type;
        mRecovery   = recovery;

	numDates = startDates.size();

        // Needs at least one protection period
        //
        if (numDates < 1)
            throw KFailure("%s: needs at least one protection period "
                           "(size = %d).\n",
                           routine,
                           numDates);

	// Check length consistency 
	ASSERT_OR_THROW(endDates.size()    == numDates);
	ASSERT_OR_THROW(settleDates.size() == numDates);

	for (idx=0; idx<=numDates-1; idx++) {
		mStartDates.insert(mStartDates.end(), startDates[idx]);
		mEndDates.insert(mEndDates.end(), endDates[idx]);
		mSettleDates.insert(mSettleDates.end(), settleDates[idx]);
		mRebates.insert(mRebates.end(), 0e0);
	}
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// General constructor with rebate.

KVPDefProtect::KVPDefProtect(
        const char            *name,         // (I) name
        const KVector(TDate)  &startDates,   // (I) Start dates
        const KVector(TDate)  &endDates,     // (I) End dates
        const KVector(TDate)  &settleDates,  // (I) Settlement dates
        const KVector(double) &rebates,      // (I) Rebates
        DefType               type,          // (I) default protection type
        double                recovery,      // (I) Recovery rate
        const char*           discZcName)    // (I) discount curve
	: KVPInstr(name)
{
static	char	routine[] = "KVPDefProtect::KVPDefProtect";

	int	idx, numDates;

    try {

	mDiscZcName = String(discZcName);
        mDefType    = type;
        mRecovery   = recovery;

	numDates = startDates.size();

        // Needs at least one protection period
        //
        if (numDates < 1)
            throw KFailure("%s: needs at least one protection period "
                           "(size = %d).\n",
                           routine,
                           numDates);

	// Check length consistency 
	ASSERT_OR_THROW(endDates.size()    == numDates);
	ASSERT_OR_THROW(settleDates.size() == numDates);
	ASSERT_OR_THROW(rebates.size()     == numDates);

	for (idx=0; idx<=numDates-1; idx++) {
		mStartDates.insert(mStartDates.end(), startDates[idx]);
		mEndDates.insert(mEndDates.end(), endDates[idx]);
		mSettleDates.insert(mSettleDates.end(), settleDates[idx]);
		mRebates.insert(mRebates.end(), rebates[idx]);
	}
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}




//---------------------------------------------------------------
void
KVPDefProtect::CheckIsValid() const
{
static	char	routine[] = "KVPDefProtect::CheckIsValid";

	int	idx, numDates;

    try {
	numDates = mStartDates.size();

        // Needs at least one protection period
        //
        if (numDates < 1)
            throw KFailure("%s: needs at least one protection period "
                           "(size = %d).\n",
                           routine,
                           numDates);

        if (mRecovery > 1e0 ||
            mRecovery < 0e0)
        {
            throw KFailure("%s: recovery rate (%lf) is not in the range "
                           "of [0, 1].\n",
                           routine,
                           mRecovery);
        }

	ASSERT_OR_THROW(mEndDates.size()    == numDates);
	ASSERT_OR_THROW(mSettleDates.size() == numDates);
	ASSERT_OR_THROW(mRebates.size()     == numDates);

	for (idx=0; idx<=numDates-1; idx++) {
            if (mStartDates[idx] > mEndDates[idx])
                throw KFailure("%s: StartDate %d (%s) > EndDate %d (%s).\n", 
                                routine,
                                idx+1, DrlTDatePrint(NULL, mStartDates[idx]),
                                idx+1, DrlTDatePrint(NULL, mEndDates[idx]));

            if (mEndDates[idx] > mSettleDates[idx])
                throw KFailure("%s: EndDate %d (%s) > SettleDate %d (%s).\n", 
                                routine,
                                idx+1, DrlTDatePrint(NULL, mEndDates[idx]),
                                idx+1, DrlTDatePrint(NULL, mSettleDates[idx]));
        }

   }
   catch (KFailure) {
        throw KFailure("%s: failed for %s.\n", routine, GetName());
   }
}




//---------------------------------------------------------------

istream&
KVPDefProtect::Get(istream& is, int drw)
{

    char   type;

    try {
	if (drw) {
	    int	idx, numDates;

	    numDates = getInt(is, "KVPDefProtect::Get: number of dates.");

	    for (idx=0; idx<numDates; idx++) {
		mStartDates.insert(mStartDates.end(),
			getTDate(is, "KVPDefProtect::Get: mStartDates."));

		mEndDates.insert(mEndDates.end(),
			getTDate(is, "KVPDefProtect::Get: mEndDates."));

		mSettleDates.insert(mSettleDates.end(),
			getTDate(is, "KVPDefProtect::Get: mSettleDates."));

		mRebates.insert(mRebates.end(),
			getDouble(is, "KVPDefProtect::Get: mRebates."));
	    }

            type = getChar(is, "KVPDefProtect::Get: mDefType.");
            if (toupper(type) == 'K')
                mDefType = DEF_KNOCKIN;
            else
                mDefType = DEF_EXPOSURE;

            mRecovery = getDouble(is, "KVPDefProtect::Get: mRecovery.");

	} else {
		throw KFailure("KVPDefProtect::Get: format N/A.\n");
	}


	return(is);
    }
    catch (KFailure) {
	throw KFailure("KVPDefProtect::Get: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPDefProtect::Put(ostream& os, int indent) const
{
	int	idx, numDates;

    try {

        CheckIsValid();

	numDates = mStartDates.size();

	os << "NAME: `" << GetName() << "'" << endl;

	os << "NUMDATES: " << numDates << endl;

	os << "    START       END       SETTLE     REBATE " << endl;
	for (idx=0; idx<numDates; idx++) {
		os << format(" %10s %10s %10s %18.6f\n",
			DrlTDatePrint(NULL, mStartDates[idx]),
			DrlTDatePrint(NULL, mEndDates[idx]),
			DrlTDatePrint(NULL, mSettleDates[idx]),
			mRebates[idx]);
	}
        if (mDefType == DEF_KNOCKIN)
	    os << "DEFAULT PROTECTION TYPE: DEF_KNOCKIN" << endl;
        else
	    os << "DEFAULT PROTECTION TYPE: DEF_EXPOSURE" << endl;

	os << "Recovery Rate: "  << mRecovery << endl;

	os << "DISCOUNT CURVE: " << mDiscZcName << endl;

	// Print dependencies
	this->KVPAtom::Put(os);

	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPDefProtect::Get: failed.\n");
    }
}



//---------------------------------------------------------------
// Only one underlying is allowed.
//
ostream&
KVPDefProtect::YacctionWrite(ostream& os, int indent)
{
	int	idx, numDates;

    try {

        CheckIsValid();

	if (GetWriteFlag())
	{

	    if (NumDep() > 1)
		throw KFailure("KVPDefProtect::YacctionWrite: Only "
			"one underlying is allowed in the wrapper.\n");

	    numDates = mStartDates.size();

            if (mDefType == DEF_KNOCKIN)
            {
	        os << GetName() << "=DEFKNOCKIN("  << endl; 
	        os << "{";

	        for (idx=0; idx<numDates; idx++) {
	            os << format(" %10s %10s %10s %15.3f\n",
		    	    GtoFormatDate(mStartDates[idx]),
			    GtoFormatDate(mEndDates[idx]),
			    GtoFormatDate(mSettleDates[idx]),
		            mRebates[idx]);
	        }
	        os << "}," << endl;
            }
            else
            {
	        os << GetName() << "=DEFEXPOSURE(" << endl; 
	        os << "{";

	        for (idx=0; idx<numDates; idx++) {
	            os << format(" %10s %10s %10s\n",
		    	    GtoFormatDate(mStartDates[idx]),
			    GtoFormatDate(mEndDates[idx]),
			    GtoFormatDate(mSettleDates[idx]));
	        }

	        os << "}," << mRecovery << "," << endl;
            }

	    os << "\t"   << Dep(0)->GetName() << "," << endl
	       << "\t\"" << mDiscZcName << "\");"
	       << endl << endl;

	    WriteDone();
	}
 
	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPDefProtect::YacctionWrite: failed.\n");
    }
}

