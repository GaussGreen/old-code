/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      D. Liu
 ************************************************************************/

#include "vpkio.h"

#include "kutilios.h"		// ios utilities

extern	"C" {
#include "datelist.h"           // GtoNewDateList 
#include "drltime.h"		// DrlTDateIntervalNil() 
#include "drlinter.h"		// DrlTDateLinearInterp1d() 
};



//-------------------------------------------------------------
// Convenience constructor to build an observation table
// from start date to and date with a given frequency.

KVPKnockIO::KVPKnockIO(
	const char *name,			// (I) name
	const KKnockIO &ioType,			// (I) Knock in/out type
	const KKnockIO &ioWindow,		// (I) Knock in/out window
	const SharedPointer<KRate> rateIndex,		// (I) Knock in/out rate index
	const KSmooth &smooth,			// (I) Node smoothing flag
	TDate startDate,			// (I) This date isn't included
	TDate matDate,				// (I) maturity date
	const TDateInterval &freq,		// (I) obs freq
	TBoolean stubAtEnd,			// (I) stub
	const KDateInterval &notDays, 		// (I) # of notif days
	
	const KVector(TDate)  &barrierDates,	// (I) barrier schedule dates
	const KVector(double) &barrierLo,	// (I) Low barries
	const KVector(double) &barrierHi,	// (I) High barries
	const KVector(double) &rebates,		// (I) Rebates
	const char*  discZcName)		// (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPKnockIO::KVPKnockIO";

	int		idx, numDates;
	TDateList	*dl = NULL;
	TDate		obsDate, settleDate;
	double		barHi, barLo, rbate;

	KRateReset	*rateReset = NULL;

	SharedPointer<KVPAtom>	koRate;

			// We need this because ALIB routines
			// do not take const 
	KDateInterval	freqA(freq);
	KVector(TDate)	barrierDatesA = barrierDates;
	KVector(double)	barrierLoA = barrierLo;
	KVector(double)	barrierHiA = barrierHi;
	KVector(double)	rebatesA = rebates;

    try {

	mObsType    = FREQ;

	mIOType     = ioType;
	mIOWindow   = ioWindow;

	mSmooth    = smooth;

	mDiscZcName = String(discZcName);

	ASSERT_OR_THROW(barrierDates.size() == barrierLo.size());
	ASSERT_OR_THROW(barrierDates.size() == barrierHi.size());
	ASSERT_OR_THROW(barrierDates.size() == rebates.size());

	
	ASSERT_OR_THROW((dl = GtoNewDateList(
		startDate,
		matDate,
		&((TDateInterval) freqA),
		stubAtEnd)) != NULL);
	numDates = dl->fNumItems;

	for (idx=0; idx<=numDates-1; idx++) {
		settleDate = dl->fArray[idx];
		obsDate    = settleDate - notDays; 


		// Rate Reset
		rateReset = new KRateReset(obsDate,
					   *rateIndex);


		mRateIndex.push_back(rateReset);
		rateReset = NULL;

		// Interpolate the barriers from arrays
		//
		IF_FAILED_THROW(DrlTDateLinearInterp1d(&barrierDatesA[0],
						       &barrierLoA[0],
						       barrierDatesA.size(),
						       obsDate,
						       &barLo)); 

		IF_FAILED_THROW(DrlTDateLinearInterp1d(&barrierDatesA[0],
						       &barrierHiA[0],
						       barrierDatesA.size(),
						       obsDate,
						       &barHi)); 

		IF_FAILED_THROW(DrlTDateLinearInterp1d(&barrierDatesA[0],
						       &rebatesA[0],
						       barrierDatesA.size(),
						       obsDate,
						       &rbate)); 

		// Quick check
		//
		if (barLo > barHi * (1. - BARRIER_TOL))
			throw KFailure("%s: low barrier (%f) > high barrier "
				       "(%f) on date %s.\n",
					routine,
					barLo,
					barHi,
					GtoFormatDate(obsDate));

		mObsDates.push_back(obsDate);
		mSettleDates.push_back(settleDate);
		mBarrierLo.push_back(barLo);
		mBarrierHi.push_back(barHi);
		mRebates.push_back(rbate);

	}


	// Add rate as dependecy
	SharedPointerConvertTo(rateIndex, koRate);
	AddDep(koRate);

	GtoFreeDateList(dl);
    }
    catch (KFailure) {
	GtoFreeDateList(dl);
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// Convenience constructor to build an observation table
// from specified dates.  Information on effective observation
// dates is contained in rateIndex's spot offset.

KVPKnockIO::KVPKnockIO(
	const char *name,			// (I) name
	const KKnockIO &ioType,			// (I) Knock in/out type
	const KKnockIO &ioWindow,		// (I) Knock in/out window
	const SharedPointer<KRate> rateIndex,		// (I) Knock in/out rate index
	const KSmooth &smooth,			// (I) Node smoothing flag

	const KVector(TDate)  &obsDates,	// (I) Observation dates 
	const KVector(TDate)  &settleDates, 	// (I) Settlement dates
	const KVector(double) &barrierLos,	// (I) Low barries
	const KVector(double) &barrierHis,	// (I) High barries
	const KVector(double) &rebates,		// (I) Rebates

	const char* 	discZcName)		// (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPKnockIO::KVPKnockIO";

    try {

	int	i;

	KRateReset	*rateReset = NULL;

	SharedPointer<KVPAtom>	koRate;

	mObsType    = DATES;

	mIOType     = ioType;
	mIOWindow   = ioWindow;

	mSmooth     = smooth;

	mDiscZcName = String(discZcName);

	ASSERT_OR_THROW(obsDates.size() == settleDates.size());
	ASSERT_OR_THROW(obsDates.size() == barrierLos.size());
	ASSERT_OR_THROW(obsDates.size() == barrierHis.size());
	ASSERT_OR_THROW(obsDates.size() == rebates.size());

	for (i=0; i<=obsDates.size()-1; i++)
	{
		if(obsDates[i] > settleDates[i])
			throw KFailure("%s: observation date (%s) > "
				       "settlement date (%s).\n",
					routine,
					GtoFormatDate(obsDates[i]),
					GtoFormatDate(settleDates[i]));

		if(barrierLos[i] > barrierHis[i]) 
			throw KFailure("%s: lower barrier value (%f) > "
				       "higher barrier value (%f).\n",
					routine,
					barrierLos[i],
					barrierHis[i]);

		// Rate Reset
		rateReset = new KRateReset(obsDates[i],
					   *rateIndex);


		mRateIndex.push_back(rateReset);
		rateReset = NULL;
	}
	
	mObsDates    = obsDates;
	mSettleDates = settleDates;
	mBarrierLo   = barrierLos;
	mBarrierHi   = barrierHis;
	mRebates     = rebates;

	// Add rate as dependecy
	SharedPointerConvertTo(rateIndex, koRate);
	AddDep(koRate);
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// Convenience constructor to build an observation table
// from specified dates.  Effective observation dates are 
// given explicitly.

KVPKnockIO::KVPKnockIO(
	const char *name,			// (I) name
	const KKnockIO &ioType,			// (I) Knock in/out type
	const KKnockIO &ioWindow,		// (I) Knock in/out window
	const SharedPointer<KRate> rateIndex,		// (I) Knock in/out rate index
	const KSmooth &smooth,			// (I) Node smoothing flag

	const KVector(TDate)  &obsDates,	// (I) Observation dates 
	const KVector(TDate)  &obsEffDates,	// (I) Observation eff dates 
	const KVector(TDate)  &settleDates, 	// (I) Settlement dates
	const KVector(double) &barrierLos,	// (I) Low barries
	const KVector(double) &barrierHis,	// (I) High barries
	const KVector(double) &rebates,		// (I) Rebates

	const char* 	discZcName)		// (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPKnockIO::KVPKnockIO";

    try {

	int	i;

	KRateReset	*rateReset = NULL;

	SharedPointer<KVPAtom>	koRate;


	mObsType    = DATES;

	mIOType     = ioType;
	mIOWindow   = ioWindow;

	mSmooth     = smooth;

	mDiscZcName = String(discZcName);

	ASSERT_OR_THROW(obsDates.size() == obsEffDates.size());
	ASSERT_OR_THROW(obsDates.size() == settleDates.size());
	ASSERT_OR_THROW(obsDates.size() == barrierLos.size());
	ASSERT_OR_THROW(obsDates.size() == barrierHis.size());
	ASSERT_OR_THROW(obsDates.size() == rebates.size());

	for (i=0; i<=obsDates.size()-1; i++)
	{
		if(obsDates[i] > settleDates[i])
			throw KFailure("%s: observation date (%s) > "
				       "settlement date (%s).\n",
					routine,
					GtoFormatDate(obsDates[i]),
					GtoFormatDate(settleDates[i]));

		if(barrierLos[i] > barrierHis[i]) 
			throw KFailure("%s: lower barrier value (%f) > "
				       "higher barrier value (%f).\n",
					routine,
					barrierLos[i],
					barrierHis[i]);

		// Rate Reset
		rateReset = new KRateReset(obsDates[i],
					   obsEffDates[i],
					   *rateIndex);


		mRateIndex.push_back(rateReset);
		rateReset = NULL;
	}
	
	mObsDates    = obsDates;
	mSettleDates = settleDates;
	mBarrierLo   = barrierLos;
	mBarrierHi   = barrierHis;
	mRebates     = rebates;

	// Add rate as dependecy
	SharedPointerConvertTo(rateIndex, koRate);
	AddDep(koRate);
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}




//-------------------------------------------------------------
// Convenience constructor to build an observation table
// from specified rates and dates.  This is the most
// general one.

KVPKnockIO::KVPKnockIO(
	const char *name,			// (I) name
	const KKnockIO &ioType,			// (I) Knock in/out type
	const KKnockIO &ioWindow,		// (I) Knock in/out window
	const KSmooth &smooth,			// (I) Node smoothing flag

	const KVector(TDate)  &obsDates,	// (I) Observation dates 
	const KVector(TDate)  &obsEffDates,	// (I) Observation dates 
	const KVector(TDate)  &settleDates, 	// (I) Settlement dates
	const KVector(SharedPointer<KRate>) &rateIndex,	// (I) Knock in/out rate index
	const KVector(double) &barrierLos,	// (I) Low barries
	const KVector(double) &barrierHis,	// (I) High barries
	const KVector(double) &rebates,		// (I) Rebates

	const char* 	discZcName)		// (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPKnockIO::KVPKnockIO";

    try {

	int	i;

	SharedPointer<KVPAtom>	koRate;

	mObsType    = DATES;

	mIOType     = ioType;
	mIOWindow   = ioWindow;

	mSmooth     = smooth;

	mDiscZcName = String(discZcName);

	ASSERT_OR_THROW(obsDates.size() == obsEffDates.size());
	ASSERT_OR_THROW(obsDates.size() == rateIndex.size());
	ASSERT_OR_THROW(obsDates.size() == settleDates.size());
	ASSERT_OR_THROW(obsDates.size() == barrierLos.size());
	ASSERT_OR_THROW(obsDates.size() == barrierHis.size());
	ASSERT_OR_THROW(obsDates.size() == rebates.size());

	for (i=0; i<=obsDates.size()-1; i++)
	{
		if(obsDates[i] > settleDates[i])
			throw KFailure("%s: observation date (%s) > "
				       "settlement date (%s).\n",
					routine,
					GtoFormatDate(obsDates[i]),
					GtoFormatDate(settleDates[i]));

		if(barrierLos[i] > barrierHis[i]) 
			throw KFailure("%s: lower barrier value (%f) > "
				       "higher barrier value (%f).\n",
					routine,
					barrierLos[i],
					barrierHis[i]);

		// Insert rate array
		//
		mRateIndex.insert(mRateIndex.end(),
					new KRateReset(
						obsDates[i],
						obsEffDates[i],
						*rateIndex[i]));

		// Add rate as dependecy
		SharedPointerConvertTo(rateIndex[i], koRate);
		AddDep(koRate);
	}
	
	mObsDates    = obsDates;
	mSettleDates = settleDates;
	mBarrierLo   = barrierLos;
	mBarrierHi   = barrierHis;
	mRebates     = rebates;
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//-------------------------------------------------------------
// Convenience constructor for tree NODES or CONTINUOUS observation types.
// 
//                          INCOMPLETE!!
//

KVPKnockIO::KVPKnockIO(
	const char 	*name,		// (I) name
	const KKnockIO	&ioType,	// (I) Knock in/out type
	const KKnockIO	&ioWindow,	// (I) Knock in/out window type
	const SharedPointer<KRate> rateIndex,	// (I) Knock in/out rate index
	const KSmooth		&smooth,	// (I) Node smoothing flag

	const KObsType	obsType,	// (I) Observation type
	const KDateInterval 	&notDays, 	// (I) # of notifcation days
	
	const KVector(TDate)	&barrierDates,	// (I) barrier schedule dates
	const KVector(double)	&barrierLo,	// (I) Low barries
	const KVector(double)	&barrierHi,	// (I) High barries
	const KVector(double)	&rebates,	// (I) Rebates

	const char* 	discZcName)     // (I) discount curve		
	: KVPInstr(name)
{
static	char	routine[] = "KVPKnockIO::KVPKnockIO";

	int		idx, numDates;
	TDate		obsDate;

	KRateReset	*rateReset = NULL;

	SharedPointer<KVPAtom>	koRate;

try {

	if (obsType == NODES || obsType == CONTINUOUS)
		mObsType    = obsType;
	else
		throw KFailure("%s: invalid observation type (%d). "
				"only NODES or CONTINUOUS is allowed.\n",
				routine,
				obsType);


	mIOType     = ioType;
	mIOWindow   = ioWindow;

	mSmooth   = smooth;

	mDiscZcName = String(discZcName);

	ASSERT_OR_THROW(barrierDates.size() == barrierLo.size());
	ASSERT_OR_THROW(barrierDates.size() == barrierHi.size());
	ASSERT_OR_THROW(barrierDates.size() == rebates.size());

	
	// NOT done properly.  Need to construct observation
	// dates based on observation type.
	//

	numDates = barrierDates.size();
	for (idx=0; idx<=numDates-1; idx++) {
		obsDate = barrierDates[idx]-notDays;
		mObsDates.push_back(obsDate);

		// Rate Reset
		rateReset = new KRateReset(obsDate,
					   *rateIndex);

		mRateIndex.push_back(rateReset);
		rateReset = NULL;

		mSettleDates.push_back(barrierDates[idx]);
		mBarrierLo.push_back(barrierLo[idx]);
		mBarrierHi.push_back(barrierHi[idx]);
		mRebates.push_back(rebates[idx]);

	}

	// Add rate as dependecy
	SharedPointerConvertTo(rateIndex, koRate);
	AddDep(koRate);
 

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}




//---------------------------------------------------------------

KVPKnockIO::~KVPKnockIO()
{
	for(KVector(KRateReset*)::iterator it=mRateIndex.begin();
			it!=mRateIndex.end(); ++it)
		delete (*it);

	mRateIndex.clear();
}



//-------------------------------------------------------------
// Include only the events with observation date >= today.

void
KVPKnockIO::ValidEvents(TDate today)		// (I) today's date
{
static	char	routine[] = "KVPKnockIO::ValidEvents";

    try {
	int	idx;

	if (mObsDates.size() <= 0) {
		throw KFailure("%s: no knock-out obs dates.\n", routine);
	}
 
	KVector(TDate)  obsDates;
	KVector(TDate)  settleDates;
	KVector(KRateReset*) rateIndex;
	KVector(double) barrierLos;
	KVector(double) barrierHis;
	KVector(double) rebates;

	for (idx=0; idx<=mObsDates.size()-1; idx++)
	{
	    if(mObsDates[idx] >= today)
            {
		obsDates.push_back(mObsDates[idx]);
		settleDates.push_back(mSettleDates[idx]);
		rateIndex.push_back(mRateIndex[idx]);
		barrierLos.push_back(mBarrierLo[idx]);
		barrierHis.push_back(mBarrierHi[idx]);
		rebates.push_back(mRebates[idx]);
            }
	    else	// free obsolate memory
		delete mRateIndex[idx];
	}
	
	// Clear up
	mObsDates.clear();
	mSettleDates.clear();
	mRateIndex.clear();
	mBarrierLo.clear();
	mBarrierHi.clear();
	mRebates.clear();

	// Assign the new list
	mObsDates    = obsDates;
	mSettleDates = settleDates;
	mRateIndex   = rateIndex;
	mBarrierLo   = barrierLos;
	mBarrierHi   = barrierHis;
	mRebates     = rebates;
	
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------

istream&
KVPKnockIO::Get(istream& is, int drw)
{

    try {
	if (drw) {
	    int	idx, numDates;
	    char smooth;

	    mIOType   = getKnockIO(is, "KVPKnockIO::Get: mIOType");
	    mIOWindow = getKnockIO(is, "KVPKnockIO::Get: mIOWindow");

	    smooth  = getChar(is, "KVPKnockIO::Get: Smooth");
	    if (toupper(smooth) == 'D')
	    	mSmooth = DOUBLE_SMOOTH;
	    else if (toupper(smooth) == 'S')
		mSmooth = SINGLE_SMOOTH;
	    else if (toupper(smooth) == 'N')
		mSmooth = NO_SMOOTH;
	    else
		throw KFailure("invalid input (%c) for smoothing.\n", 
				smooth);

	    mObsType = getKObsType(is, "KVPKnockIO::Get: mObsType");

	    numDates = getInt(is, "KVPKnockIO::Get: number of dates.");

	    for (idx=0; idx<numDates; idx++) {
		mObsDates.push_back(
			getTDate(is, "KVPKnockIO::Get: mObsDates."));

		mSettleDates.push_back(
			getTDate(is, "KVPKnockIO::Get: mSettleDates."));

		mBarrierLo.push_back(
			getDouble(is, "KVPKnockIO::Get: mBarrierLo."));

		mBarrierHi.push_back(
			getDouble(is, "KVPKnockIO::Get: mBarrierHi."));

		mRebates.push_back(
			getDouble(is, "KVPKnockIO::Get: mRebate."));
	    }

	} else {
		throw KFailure("KVPKnockIO::Get: format N/A.\n");
	}


	return(is);
    }
    catch (KFailure) {
	throw KFailure("KVPKnockIO::Get: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPKnockIO::Put(ostream& os, int indent) const
{
	int	idx, numDates;

    try {
	numDates = mObsDates.size();
	ASSERT_OR_THROW(mSettleDates.size() == numDates);
	ASSERT_OR_THROW(mBarrierLo.size() == numDates);
	ASSERT_OR_THROW(mBarrierHi.size() == numDates);
	ASSERT_OR_THROW(mRebates.size() == numDates);

	os << "NAME: `" << GetName() << "'" << endl;
	os << "KNOCK IN/OUT: "  << mIOType   << endl;
	os << "KNOCK WINDOW : " << mIOWindow << endl;

	if (mSmooth == DOUBLE_SMOOTH)
		os << "SMOOTHING: DOUBLE" << endl;
	else if (mSmooth == SINGLE_SMOOTH)
		os << "SMOOTHING: SINGLE" << endl;
	else if (mSmooth == NO_SMOOTH)
		os << "SMOOTHING: NO" << endl;
	else
		throw KFailure("invalid smoothing method (%d).\n", mSmooth);
			

	os << "NUMDATES: " << numDates << endl;

	os << "  INDEX   OBSERVE      SETTLE        BARRIER_LO   BARRIER_HI   REBATE" << endl;
	for (idx=0; idx<numDates; idx++) {
		os << mRateIndex[idx]->Rate().Maturity()
		   << " " << mRateIndex[idx]->Rate().CurveName();
		
		os << format(" %10s %10s %12.6f %12.6f %12.6f\n",
			DrlTDatePrint(NULL, mObsDates[idx]),
			DrlTDatePrint(NULL, mSettleDates[idx]),
			mBarrierLo[idx],
			mBarrierHi[idx],
			mRebates[idx]);
	}
	os << "DISCOUNT CURVE: " << mDiscZcName << endl;

	// Print dependencies
	this->KVPAtom::Put(os);

	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPKnockIO::Get: failed.\n");
    }
}



//---------------------------------------------------------------
//
/*
void
KVPKnockIO::SetWriteFlag(bool flag)
{
	int     idx, numRates;
 
	KVPAtom::SetWriteFlag(flag);
 
	// Set the rate flag
	numRates = mResetIndex.size();
        
	for (idx=0; idx<numRates; idx++) {
		mResetIndex[idx]->Rate().SetWriteFlag(flag);    
	}
               
}
*/



//---------------------------------------------------------------

ostream&
KVPKnockIO::YacctionWrite(ostream& os, int indent)
{
	int	idx, numDates;

	bool	isOneRate = true;
    try {
	
	if (GetWriteFlag())
	{
	    numDates = mObsDates.size();
	    ASSERT_OR_THROW(mSettleDates.size() == numDates);
	    ASSERT_OR_THROW(mBarrierLo.size() == numDates);
	    ASSERT_OR_THROW(mBarrierHi.size() == numDates);
	    ASSERT_OR_THROW(mRebates.size() == numDates);


	    // Write the floating rate first
	    // Test if all the rates are the same for each reset
	    //
	    for (idx=1; idx<numDates; idx++) {
		if (!(mRateIndex[0]->Rate() == mRateIndex[idx]->Rate()))
		isOneRate = false;
	    }

		
	    if (!isOneRate)
		throw KFailure("KVPKnockIO::YacctionWrite: different trigger "
			"rates are not allowed in the wrapper.\n");

	    os << GetName() << "=KNOCKIO("  << endl
	       << "\t\"" << (mIOType==CRX_KNOCK_IN ? "I" : "O") << "\"," 
	       << "\"" << (mIOWindow==CRX_KNOCK_IN ? "I" : "O") << "\"," 
	       << endl;

	    if (isOneRate)
	    {
		for (idx=0; idx <NumDep(); idx++)
		{
		    if (Dep(idx)->IsType("KRate"))
		    {
	    	    	os << "\t" << Dep(idx)->GetName() << ", " << endl;
			break;
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
		   << formatDouble(" %12.6f", mRebates[idx]);

		if (!isOneRate)
		{                        
			os << " ";
			mRateIndex[idx]->Rate().WriteSimple(os);
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
	throw KFailure("KVPKnockIO::YacctionWrite: failed.\n");
    }
}

