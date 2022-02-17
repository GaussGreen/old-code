/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      D. Liu
 ************************************************************************/

#include "vpkio2idx.h"

#include "kutilios.h"		// ios utilities

extern	"C" {
#include "datelist.h"           // GtoNewDateList 
#include "drltime.h"		// DrlTDateIntervalNil() 
#include "drlinter.h"		// DrlTDateLinearInterp1d() 
};


//-------------------------------------------------------------
// Convenience constructor given observation freqency.
// Date-independent dual ko rate index.

KVPKnockIO2Idx::KVPKnockIO2Idx(
    const char  *name,                    // (I) Name
    const KKnockIO  &ioType,              // (I) Knock in/out type
    const KKnockIO  &ioWindow1,           // (I) Pay in/out 1
    const SharedPointer<KRate> rateIndex1,// (I) Knock in/out rate index 1
    const KKnockIO  &ioWindow2,           // (I) Pay in/out 2
    const SharedPointer<KRate> rateIndex2,// (I) Knock in/out rate index 2
    const KSmooth   &smooth,              // (I) Node smoothing flag
    TDate           startDate,            // (I) Not included
    TDate           matDate,              // (I) Maturity date
    const TDateInterval &freq,            // (I) Obs frequency as interval
    TBoolean        stubAtEnd,            // (I) stub
    const KDateInterval   &notDays,       // (I) # of notif days
    const KVector(TDate)  &barrierDates,  // (I) Barrier dates
    const KVector(double) &barrierLo1,    // (I) Low barrier 1
    const KVector(double) &barrierHi1,    // (I) High barrier 1
    const KVector(double) &barrierLo2,    // (I) Low barrier 2
    const KVector(double) &barrierHi2,    // (I) High barrier 2
    const KVector(double) &rebates,       // (I) Rebates
    const char* discZcName)               // (I) discount curve
:KVPKnockIO(
    name,
    ioType,
    ioWindow1,
    rateIndex1,
    smooth,
    startDate,
    matDate,
    freq,
    stubAtEnd,
    notDays,
    barrierDates,
    barrierLo1,
    barrierHi1,
    rebates,
    discZcName)
{
static  char    routine[] = "KVPKnockIO2Idx::KVPKnockIO2Idx";

    int         idx;
    TDate       obsDate;
    double      barHi2, barLo2;

    KRateReset  *rateReset2 = NULL;

    SharedPointer<KVPAtom>  koRate2;

            // We need this because ALIB routines
            // do not take const 
    KVector(TDate)  barrierDatesA = barrierDates;
    KVector(double) barrierLoA = barrierLo2;
    KVector(double) barrierHiA = barrierHi2;

    try {

    mIOWindow2   = ioWindow2;

    ASSERT_OR_THROW(barrierDates.size() == barrierLo2.size());
    ASSERT_OR_THROW(barrierDates.size() == barrierHi2.size());

    for (idx=0; idx<=mObsDates.size()-1; idx++) {
        obsDate = mObsDates[idx];

        // Rate reset -- the reset effective is embedded 
        // in the rate definition
        rateReset2 = new KRateReset(obsDate,
                                    *rateIndex2);

        mRateIndex2.push_back(rateReset2);
        rateReset2 = NULL;

        // Interpolate the barriers from arrays
        //
        IF_FAILED_THROW(DrlTDateLinearInterp1d(&barrierDatesA[0],
                               &barrierLoA[0],
                               barrierDatesA.size(),
                               obsDate,
                               &barLo2)); 

        IF_FAILED_THROW(DrlTDateLinearInterp1d(&barrierDatesA[0],
                               &barrierHiA[0],
                               barrierDatesA.size(),
                               obsDate,
                               &barHi2)); 

        // Quick check
        //
        if (barLo2 > barHi2 * (1. - BARRIER_TOL))
            throw KFailure("%s: low barrier 2 (%f) > "
                           "high barrier 2 (%f) on date %s.\n",
                    routine,
                    barLo2,
                    barHi2,
                    GtoFormatDate(obsDate));

        mObsDates2.push_back(obsDate);
        mBarrierLo2.push_back(barLo2);
        mBarrierHi2.push_back(barHi2);
    }


    // Add rate as dependecy
    SharedPointerConvertTo(rateIndex2, koRate2);
    AddDep(koRate2);

    }
    catch (KFailure) {
        throw KFailure("%s: failed for %s.\n", routine, name);
    }
}


//-------------------------------------------------------------
// Convenience constructor for arbitrary observation and
// settlement dates.  Information on effective observation dates 
// implied from the rateIndex spot offset. 
// Time-independent dual ko rate index.
KVPKnockIO2Idx::KVPKnockIO2Idx(
	const char *name,		    // (I) name
	const KKnockIO &ioType,		// (I) Knock in/out type
	const KKnockIO &ioWindow,	// (I) Pay in/out
	const SharedPointer<KRate> &rateIndex, // (I) Ko rate index
	const KKnockIO &ioWindow2,	// (I) Pay in/out
	const SharedPointer<KRate> &rateIndex2,// (I) Ko rate index
	const KSmooth &smooth,      // (I) Node smoothing flag

	const KVector(TDate)  &obsDates,      // (I) Observ dates
	const KVector(TDate)  &settleDates,   // (I) Settlement dates
	const KVector(double) &barrierLo,     // (I) Low barrier  1
	const KVector(double) &barrierHi,     // (I) High barrier 1
	const KVector(double) &barrierLo2,    // (I) Low barrier  2
	const KVector(double) &barrierHi2,    // (I) High barrier 2
	const KVector(double) &rebates,       // (I) Rebates
	const char*     discZcName)           // (I) discount curve
: KVPKnockIO(
    name,
    ioType,
    ioWindow,
    rateIndex,
    smooth,
    obsDates,
    settleDates,
    barrierLo,
    barrierHi,
    rebates,
    discZcName)
{
static  char    routine[] = "KVPKnockIO2Idx::KVPKnockIO2Idx";

    int         idx;
    TDate       obsDate;

    KRateReset  *rateReset2 = NULL;

    SharedPointer<KVPAtom>  koRate2;

    try {

    mIOWindow2   = ioWindow2;

    ASSERT_OR_THROW(obsDates.size() == barrierLo2.size());
    ASSERT_OR_THROW(obsDates.size() == barrierHi2.size());

    for (idx=0; idx<=mObsDates.size()-1; idx++) {
        obsDate    = mObsDates[idx];

        // Rate reset -- the reset effective is input directly 
        rateReset2 = new KRateReset(obsDate,
                                    *rateIndex2);

        mRateIndex2.push_back(rateReset2);
        rateReset2 = NULL;

        // Create duplicate observation date list 
        mObsDates2.push_back(obsDate);
    }

    mBarrierLo2 = barrierLo2;
    mBarrierHi2 = barrierHi2;

    // Add rate as dependecy
    SharedPointerConvertTo(rateIndex2, koRate2);
    AddDep(koRate2);

    }
    catch (KFailure) {
        throw KFailure("%s: failed for %s.\n", routine, name);
    }
}


//-------------------------------------------------------------
// Convenience constructor for arbitrary observation and
// settlement dates.  Information on effective observation dates 
// is given explicitly. Time-independent dual ko rate index.
//
KVPKnockIO2Idx::KVPKnockIO2Idx(
	const char *name,		    // (I) name
	const KKnockIO &ioType,		// (I) Knock in/out type
	const KKnockIO &ioWindow,	// (I) Pay in/out
	const SharedPointer<KRate> &rateIndex, // (I) Ko rate index
	const KKnockIO &ioWindow2,	// (I) Pay in/out
	const SharedPointer<KRate> &rateIndex2,// (I) Ko rate index
	const KSmooth &smooth,      // (I) Node smoothing flag

	const KVector(TDate)  &obsDates,      // (I) Observ dates
	const KVector(TDate)  &obsEffDates,   // (I) Observ eff dates
	const KVector(TDate)  &settleDates,   // (I) Settlement dates
	const KVector(double) &barrierLo,     // (I) Low barrier  1
	const KVector(double) &barrierHi,     // (I) High barrier 1
	const KVector(double) &barrierLo2,    // (I) Low barrier  2
	const KVector(double) &barrierHi2,    // (I) High barrier 2
	const KVector(double) &rebates,       // (I) Rebates
	const char*     discZcName)           // (I) discount curve
: KVPKnockIO(
    name,
    ioType,
    ioWindow,
    rateIndex,
    smooth,
    obsDates,
    obsEffDates,
    settleDates,
    barrierLo,
    barrierHi,
    rebates,
    discZcName)
{
static  char    routine[] = "KVPKnockIO2Idx::KVPKnockIO2Idx";

    int         idx;
    TDate       obsDate, obsEffDate;

    KRateReset  *rateReset2 = NULL;

    SharedPointer<KVPAtom>  koRate2;

    try {

    mIOWindow2   = ioWindow2;

    ASSERT_OR_THROW(obsDates.size() == barrierLo2.size());
    ASSERT_OR_THROW(obsDates.size() == barrierHi2.size());

    for (idx=0; idx<=mObsDates.size()-1; idx++) {
        obsDate    = mObsDates[idx];
        obsEffDate = obsEffDates[idx];

        // Rate reset -- the reset effective is input directly 
        rateReset2 = new KRateReset(obsDate,
                                    obsEffDate,
                                    *rateIndex2);

        mRateIndex2.push_back(rateReset2);
        rateReset2 = NULL;

        mObsDates2.push_back(obsDate);
    }

    mBarrierLo2 = barrierLo2;
    mBarrierHi2 = barrierHi2;

    // Add rate as dependecy
    SharedPointerConvertTo(rateIndex2, koRate2);
    AddDep(koRate2);

    }
    catch (KFailure) {
        throw KFailure("%s: failed for %s.\n", routine, name);
    }
}



//-------------------------------------------------------------
// Convenience constructor for arbitrary observation dates
// settlement dates.  Effective observation dates are given
// explicitly.

KVPKnockIO2Idx::KVPKnockIO2Idx(
    const char *name,           // (I) name
    const KKnockIO &ioType,     // (I) Knock in/out type
    const KKnockIO &ioWindow,   // (I) Pay in/out
    const KVector(SharedPointer<KRate>) &rateIndex,//  (I) Ko rate index
    const KKnockIO &ioWindow2,  // (I) Pay in/out
    const KVector(SharedPointer<KRate>) &rateIndex2,// (I) Ko rate index
    const KSmooth &smooth,      // (I) Node smoothing flag
    const KVector(TDate)  &obsDates,      // (I) Observ dates
    const KVector(TDate)  &obsEffDates,   // (I) Observ eff dates
    const KVector(TDate)  &settleDates,   // (I) Settlement dates
    const KVector(double) &barrierLo,     // (I) Low barrier  1
    const KVector(double) &barrierHi,     // (I) High barrier 1
    const KVector(double) &barrierLo2,    // (I) Low barrier  2
    const KVector(double) &barrierHi2,    // (I) High barrier 2
    const KVector(double) &rebates,       // (I) Rebates
    const char*     discZcName)           // (I) discount curve
: KVPKnockIO(
    name,
    ioType,
    ioWindow,
    smooth,
    obsDates,
    obsEffDates,
    settleDates,
    rateIndex,
    barrierLo,
    barrierHi,
    rebates,
    discZcName)
{
static  char    routine[] = "KVPKnockIO2Idx::KVPKnockIO2Idx";

    int         i;
    TDate       obsDate, obsEffDate;

    KRateReset  *rateReset2 = NULL;

    SharedPointer<KVPAtom>  koRate2;

    try {

    mIOWindow2   = ioWindow2;

    ASSERT_OR_THROW(obsDates.size() == barrierLo2.size());
    ASSERT_OR_THROW(obsDates.size() == barrierHi2.size());
    ASSERT_OR_THROW(obsDates.size() == rateIndex2.size());

    for (i=0; i<=mObsDates.size()-1; i++) 
    {
        obsDate    = mObsDates[i];
        obsEffDate = obsEffDates[i];

        // Rate reset -- the reset effective is input directly 
        rateReset2 = new KRateReset(obsDate,
                                    obsEffDate,
                                    *rateIndex2[i]);

        mRateIndex2.insert(mRateIndex2.end(), rateReset2);
        rateReset2 = NULL;

        // Add rate as dependecy
        SharedPointerConvertTo(rateIndex2[i], koRate2);
        AddDep(koRate2);

        mObsDates2.push_back(obsDate);
    }

    mBarrierLo2 = barrierLo2;
    mBarrierHi2 = barrierHi2;

    }
    catch (KFailure) {
        throw KFailure("%s: failed for %s.\n", routine, name);
    }
}






//---------------------------------------------------------------

KVPKnockIO2Idx::~KVPKnockIO2Idx()
{
    for(KVector(KRateReset*)::iterator it=mRateIndex2.begin();
            it!=mRateIndex2.end(); ++it)
        delete (*it);

    mRateIndex2.clear();
}



//-------------------------------------------------------------
// Include only the events with observation date >= today.

void
KVPKnockIO2Idx::ValidEvents(TDate today)    // (I) today's date
{
static  char    routine[] = "KVPKnockIO2Idx::ValidEvents";

    try {
    int idx;

    if (mObsDates.size() <= 0) {
        throw KFailure("%s: no knock-out obs dates.\n", routine);
    }

    KVector(TDate)  obsDates2;
    KVector(KRateReset*) rateIndex2;
    KVector(double) barrierLos2;
    KVector(double) barrierHis2;

    for (idx=0; idx<=mObsDates2.size()-1; idx++)
    {
        if(mObsDates2[idx] >= today)
        {
            obsDates2.push_back(mObsDates2[idx]);
            rateIndex2.push_back(mRateIndex2[idx]);
            barrierLos2.push_back(mBarrierLo2[idx]);
            barrierHis2.push_back(mBarrierHi2[idx]);
        }
        else // free obsolete memory
        {
            delete mRateIndex2[idx];
        }
    }
    
    // Clear up
    mObsDates2.clear();
    mRateIndex2.clear();
    mBarrierLo2.clear();
    mBarrierHi2.clear();

    // Assign the new list
    mObsDates2   = obsDates2;
    mRateIndex2  = rateIndex2;
    mBarrierLo2  = barrierLos2;
    mBarrierHi2  = barrierHis2;
    
    }
    catch (KFailure) {
    throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------

istream&
KVPKnockIO2Idx::Get(istream& is, int drw)
{

    try {
    if (drw) {
        int idx, numDates;
        char smooth;

        mIOType    = getKnockIO(is, "KVPKnockIO2Idx::Get: mIOType");
        mIOWindow  = getKnockIO(is, "KVPKnockIO2Idx::Get: mIOWindow1");
        mIOWindow2 = getKnockIO(is, "KVPKnockIO2Idx::Get: mIOWindow2");

        smooth  = getChar(is, "KVPKnockIO2Idx::Get: Smooth");
        if (toupper(smooth) == 'D')
            mSmooth = DOUBLE_SMOOTH;
        else if (toupper(smooth) == 'S')
        mSmooth = SINGLE_SMOOTH;
        else if (toupper(smooth) == 'N')
        mSmooth = NO_SMOOTH;
        else
        throw KFailure("invalid input (%c) for smoothing.\n", 
                smooth);

        mObsType = getKObsType(is, "KVPKnockIO2Idx::Get: mObsType");

        numDates = getInt(is, "KVPKnockIO2Idx::Get: number of dates.");

        for (idx=0; idx<numDates; idx++) {
        mObsDates.push_back(
            getTDate(is, "KVPKnockIO2Idx::Get: mObsDates."));

        mSettleDates.push_back(
            getTDate(is, "KVPKnockIO2Idx::Get: mSettleDates."));

        mBarrierLo.push_back(
            getDouble(is, "KVPKnockIO2Idx::Get: mBarrierLo."));

        mBarrierHi.push_back(
            getDouble(is, "KVPKnockIO2Idx::Get: mBarrierHi."));

       mBarrierLo2.push_back(
            getDouble(is, "KVPKnockIO2Idx::Get: mBarrierLo2."));

        mBarrierHi2.push_back(
            getDouble(is, "KVPKnockIO2Idx::Get: mBarrierHi2."));

        mRebates.push_back(
            getDouble(is, "KVPKnockIO2Idx::Get: mRebate."));
        }

    } else {
        throw KFailure("KVPKnockIO2Idx::Get: format N/A.\n");
    }


    return(is);
    }
    catch (KFailure) {
    throw KFailure("KVPKnockIO2Idx::Get: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPKnockIO2Idx::Put(ostream& os, int indent) const
{
    int idx, numDates;

    try {
    numDates = mObsDates.size();
    ASSERT_OR_THROW(mSettleDates.size() == numDates);
    ASSERT_OR_THROW(mBarrierLo.size()   == numDates);
    ASSERT_OR_THROW(mBarrierHi.size()   == numDates);
    ASSERT_OR_THROW(mBarrierLo2.size()  == numDates);
    ASSERT_OR_THROW(mBarrierHi2.size()  == numDates);
    ASSERT_OR_THROW(mRebates.size()     == numDates);

    os << "NAME: `" << GetName() << "'" << endl;
    os << "KNOCK IN/OUT: "   << mIOType   << endl;
    os << "KNOCK WINDOW1 : " << mIOWindow << endl;
    os << "KNOCK WINDOW2 : " << mIOWindow2<< endl;

    if (mSmooth == DOUBLE_SMOOTH)
        os << "SMOOTHING: DOUBLE" << endl;
    else if (mSmooth == SINGLE_SMOOTH)
        os << "SMOOTHING: SINGLE" << endl;
    else if (mSmooth == NO_SMOOTH)
        os << "SMOOTHING: NO" << endl;
    else
        throw KFailure("invalid smoothing method (%d).\n", mSmooth);
            

    os << "NUMDATES: " << numDates << endl;

    os << "  INDEX1   INDEX2   OBSERVE      SETTLE        BARRIER_LO1  "
       << "BARRIER_HI1  BARRIER_LO2  BARRIER_HI2  REBATE" << endl;
    for (idx=0; idx<numDates; idx++) {
        os << mRateIndex[idx]->Rate().Maturity()
           << " " << mRateIndex[idx]->Rate().CurveName();
        os << mRateIndex2[idx]->Rate().Maturity()
           << " " << mRateIndex2[idx]->Rate().CurveName();
        
        os << format(" %10s %10s %12.6f %12.6f %12.6f\n",
            DrlTDatePrint(NULL, mObsDates[idx]),
            DrlTDatePrint(NULL, mSettleDates[idx]),
            mBarrierLo[idx],
            mBarrierHi[idx],
            mBarrierLo2[idx],
            mBarrierHi2[idx],
            mRebates[idx]);
    }
    os << "DISCOUNT CURVE: " << mDiscZcName << endl;

    // Print dependencies
    this->KVPAtom::Put(os);

    return(os);
    }
    catch (KFailure) {
    throw KFailure("KVPKnockIO2Idx::Put: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPKnockIO2Idx::YacctionWrite(ostream& os, int indent)
{
    int idx, numDates;

    bool    isOneRate = true;
    try {
    
    if (GetWriteFlag())
    {
        numDates = mObsDates.size();
        ASSERT_OR_THROW(mSettleDates.size() == numDates);
        ASSERT_OR_THROW(mBarrierLo.size()   == numDates);
        ASSERT_OR_THROW(mBarrierHi.size()   == numDates);
        ASSERT_OR_THROW(mBarrierLo2.size()  == numDates);
        ASSERT_OR_THROW(mBarrierHi2.size()  == numDates);
        ASSERT_OR_THROW(mRebates.size()     == numDates);


        // Write the floating rate first
        // Test if all the rates are the same for each reset
        //
        for (idx=1; idx<numDates; idx++) {
        if (!(mRateIndex [0]->Rate() == mRateIndex [idx]->Rate()) ||
            !(mRateIndex2[0]->Rate() == mRateIndex2[idx]->Rate()) )
        isOneRate = false;
        }

        
        if (!isOneRate)
        throw KFailure("KVPKnockIO2Idx::YacctionWrite: different trigger "
            "rates are not allowed in the wrapper.\n");

        os << GetName() << "=KNOCKIO("  << endl
           << "\t\"" << (mIOType ==CRX_KNOCK_IN ? "I" : "O") << "\"," 
           << "\"" << (mIOWindow ==CRX_KNOCK_IN ? "I" : "O") << "\"," 
           << "\"" << (mIOWindow2==CRX_KNOCK_IN ? "I" : "O") << "\"," 
           << endl;

        if (isOneRate)
        {
        for (idx=0; idx <NumDep(); idx++)
        {
            if (Dep(idx)->IsType("KRate"))
            {
                os << "\t" << Dep(idx)->GetName() << ", " << endl;
            }
        }
        }
    
        if (mSmooth == DOUBLE_SMOOTH)
        os << "\t\"DOUBLE\", " << endl;
        else if (mSmooth == SINGLE_SMOOTH)
        os << "\t\"SINGLE\", " << endl;
        else if (mSmooth == NO_SMOOTH)
        os << "\t\"NONE\", " << endl;
        else
        throw KFailure("invalid smoothing method (%d).\n", mSmooth);
            
        os << "{"; 
        for (idx=0; idx<numDates; idx++) {
        
        os << format(" %10s %10s %10s",
            GtoFormatDate(mObsDates[idx]),
            GtoFormatDate(mRateIndex[idx]->EffDate()),
            GtoFormatDate(mSettleDates[idx]));

        os << formatDouble(" %12.6f", mBarrierLo[idx])
           << formatDouble(" %12.6f", mBarrierHi[idx])
           << formatDouble(" %12.6f", mBarrierLo2[idx])
           << formatDouble(" %12.6f", mBarrierHi2[idx])
           << formatDouble(" %12.6f", mRebates[idx]);

        if (!isOneRate)
        {                        
            os << " ";
            mRateIndex[idx]->Rate().WriteSimple(os);
            mRateIndex2[idx]->Rate().WriteSimple(os);
        }

        os << endl;
        }

        os << "}," << endl;
        
        for (idx=0; idx <NumDep(); idx++)
        {
        if (!Dep(idx)->IsType("KRate"))
            os << "\t" << Dep(idx)->GetName() << "," << endl;
        }

        os << "\t\"" << mDiscZcName << "\");"
           << endl << endl;

        WriteDone();    
    }

    return(os);
    }
    catch (KFailure) {
    throw KFailure("KVPKnockIO2Idx::YacctionWrite: failed.\n");
    }
}

