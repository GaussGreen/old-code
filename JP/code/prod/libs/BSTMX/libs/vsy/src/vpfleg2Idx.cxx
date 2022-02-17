/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:  David Liu
 ************************************************************************/
#include "vpfleg2Idx.h"

#include "kutilios.h"		// ios utilities


extern	"C" {
#include "datelist.h"           // GtoNewDateList 
#include "drltime.h"		// DrlTDatePrint

};



//---------------------------------------------------------------


KVPFloatLeg2Idx::KVPFloatLeg2Idx(
	const char *name,		// (I) name
	TDate startDate,		// (I) This date is not included
	TDate matDate,			// (I) maturity date
	const KDateInterval &freq,	// (I) frequency as interval
	KDayCc dayCc,			// (I) 
	KStubConv stubConv,		// (I) 
	TBoolean stubAtEnd,		// (I) T=stub at end, F=at beg 
	const SharedPointer<KRate> floatRate1, // (I) 1st rate index
	const SharedPointer<KRate> floatRate2, // (I) 2nd rate index
	const char 	*formula,	// (I)
	const char 	*discZcName)	// (I) 
  : KVPFloatLeg(
	name,		
	startDate,
	matDate,
	freq,
	dayCc,
	stubConv,
	stubAtEnd,
	floatRate1,
	formula,
	discZcName) 
{
static	char	routine[] = "KVPFloatLeg2Idx::KVPFloatLeg2Idx";
	int		idx;

	SharedPointer<KVPAtom>  rateIndex2;

 try {

	// Insert the 2nd rate index, assuming same reset date,
	// payment date and accural dates as the 1st rate index.
	//
	for (idx=0; idx < mResetDates.size(); idx++) {
		mRateResets2.insert(mRateResets2.end(),
				   new KRateReset(
						mResetDates[idx],
						mResetEffDates[idx],
						*floatRate2)
				   );
	}


	// Add rate as dependecy
	SharedPointerConvertTo(floatRate2, rateIndex2);
	AddDep(rateIndex2);

    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}




//---------------------------------------------------------------


KVPFloatLeg2Idx::KVPFloatLeg2Idx(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effecitve dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const SharedPointer<KRate> floatRate1, // (I) 1st rate index
	const SharedPointer<KRate> floatRate2, // (I) 2nd rate index
	const char *formula,			// (I) paymt formula (or NULL)
	const char *discZcName)			// (I) discount curve name
  : KVPFloatLeg(
	name,
	resetDates,
	resetEffDates,
	accStartDates,
	accEndDates,
	payDates,
	notionals,
	dayCc,	
	stubConv,
	floatRate1,
	formula,
	discZcName)
{
static	char	routine[] = "KVPFloatLeg2Idx::KVPFloatLeg2Idx";
	int	numDates, idx;

	SharedPointer<KVPAtom>  rateIndex2;

 try {

	numDates = resetDates.size();

	for (idx=0; idx < numDates; idx++) {
	   	 mRateResets2.insert(mRateResets2.end(),
			             new KRateReset(
						resetDates[idx],
						resetEffDates[idx],
						*floatRate2)
			      );
	}


	// Add rate as dependecy
	SharedPointerConvertTo(floatRate2, rateIndex2);
	AddDep(rateIndex2);

    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}


//---------------------------------------------------------------


KVPFloatLeg2Idx::KVPFloatLeg2Idx(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effecitve dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	const KVector(String)& formula,		// (I) payment formula array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const SharedPointer<KRate> floatRate1, // (I) 1st rate index
	const SharedPointer<KRate> floatRate2, // (I) 2nd rate index
	const char *discZcName)			// (I) discount curve name
  : KVPFloatLeg(
	name,
	resetDates,
	resetEffDates,
	accStartDates,
	accEndDates,
	payDates,
	notionals,
	formula,
	dayCc,	
	stubConv,
	floatRate1,
	discZcName)
{
static	char	routine[] = "KVPFloatLeg2Idx::KVPFloatLeg2Idx";
	int	numDates, idx;

	SharedPointer<KVPAtom>  rateIndex2;

 try {

	numDates = resetDates.size();

	for (idx=0; idx < numDates; idx++) {
	   	 mRateResets2.insert(mRateResets2.end(),
			             new KRateReset(
						resetDates[idx],
						resetEffDates[idx],
						*floatRate2)
			      );
	}


	// Add rate as dependecy
	SharedPointerConvertTo(floatRate2, rateIndex2);
	AddDep(rateIndex2);
    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}



//---------------------------------------------------------------


KVPFloatLeg2Idx::KVPFloatLeg2Idx(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effecitve dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	const KVector(String)& formula,		// (I) payment formula array
	const KVector(SharedPointer<KRate>) &floatRate1,// (I) paid rate 1
	const KVector(SharedPointer<KRate>) &floatRate2,// (I) paid rate 2
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const char *discZcName)			// (I) discount curve name
  : KVPFloatLeg(
	name,
	resetDates,
	resetEffDates,
	accStartDates,
	accEndDates,
	payDates,
	notionals,
	formula,
	floatRate1,
	dayCc,	
	stubConv,
	discZcName)
{
static	char	routine[] = "KVPFloatLeg2Idx::KVPFloatLeg2Idx";
	int	numDates, idx;

	SharedPointer<KVPAtom>  rateIndex2;

 try {

	numDates = resetDates.size();

	for (idx=0; idx < numDates; idx++) {
		// Add rate as dependecy
		SharedPointerConvertTo(floatRate2[idx], rateIndex2);
		AddDep(rateIndex2);

		mRateResets2.insert(mRateResets2.end(),
			             new KRateReset(
						resetDates[idx],
						resetEffDates[idx],
						*floatRate2[idx])
			      );
	}

    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}




//---------------------------------------------------------------
//
KVPFloatLeg2Idx::~KVPFloatLeg2Idx()
{
	for(KVector(KRateReset*)::iterator it=mRateResets2.begin();
                        it!=mRateResets2.end(); ++it)
                delete (*it);
 
        mRateResets2.clear();

}



//---------------------------------------------------------------
// Only called after construction

void
KVPFloatLeg2Idx::CheckValid() const
{
static	char	routine[] = "KVPFloatLeg2Idx::CheckValid";
    int	     idx, numDates;
    TDate    resetDate, accStartDate, accEndDate, payDate;

    try {
	numDates = mResetDates.size();

	ASSERT_OR_THROW(mResetDates.size()    == numDates);
	ASSERT_OR_THROW(mResetEffDates.size() == numDates);
	ASSERT_OR_THROW(mAccStartDates.size() == numDates);
	ASSERT_OR_THROW(mAccEndDates.size()   == numDates);
	ASSERT_OR_THROW(mPayDates.size()      == numDates);
	ASSERT_OR_THROW(mNotionals.size()     == numDates);
	ASSERT_OR_THROW(mRateResets.size()    == numDates);
	ASSERT_OR_THROW(mFormulas.size()      == numDates);

	ASSERT_OR_THROW(mRateResets2.size()   == numDates);


	for (idx=0; idx < numDates; idx++) {
            resetDate      = mResetDates[idx];
            accStartDate   = mAccStartDates[idx];
            accEndDate     = mAccEndDates[idx];
            payDate        = mPayDates[idx];

            if (accStartDate > accEndDate)
                throw KFailure("%s: accStartDate %d (%s) >"
				" accEndDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accStartDate),
				idx+1, DrlTDatePrint(NULL, accEndDate));

            if (resetDate > payDate)
                throw KFailure("%s: resetDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, resetDate),
				idx+1, DrlTDatePrint(NULL, payDate));

            if (accEndDate > payDate)
                throw KFailure("%s: accEndDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accEndDate),
				idx+1, DrlTDatePrint(NULL, payDate));


            //
            // Reset in arrear is NOT allowed for BOND stub
            //
            if ( (long)mStubConv == GTO_STUB_BOND &&
                 accStartDate < resetDate
               )
               throw KFailure("%s: reset in arrear is not allowed "
                              "for BOND stub "
                              "(reset date: %s, accural period: [%s, %s] ).\n",
                              routine,
                              DrlTDatePrint(NULL, resetDate),
                              DrlTDatePrint(NULL, accStartDate),
                              DrlTDatePrint(NULL, accEndDate));

            //
            // Reset in arrear is NOT allowed for default accural
            //
            if ( mDayCc.isNegative() &&
                 accStartDate < resetDate
               )
               throw KFailure("%s: reset in arrear is not allowed "
                              "for default accrual convention "
                              "(reset date: %s, accural period: [%s, %s] ).\n",
                              routine,
                              DrlTDatePrint(NULL, resetDate),
                              DrlTDatePrint(NULL, accStartDate),
                              DrlTDatePrint(NULL, accEndDate));

        }
    }
    catch (KFailure) {
	throw KFailure("KVPFloatLeg2Idx::CheckValid: failed for %s.\n", GetName());
    }

}


//---------------------------------------------------------------

ostream&
KVPFloatLeg2Idx::Put(ostream& os, int indent) const
{
	int	idx, numDates;

    try {
	numDates = mResetDates.size();

    CheckValid();

	os << "NAME: `" << GetName() << "'" << endl;

    /* ALIB Day Counts are positive integers */
    if (!mDayCc.isNegative())
        os << "DAYCC: " << mDayCc << endl;
    else
        os << "DAYCC: " << mDayCc << "(Accrual on default)" << endl;

	os << "STUBCONV: " << mStubConv << endl;
	os << "DISCOUNT CURVE: " << mDiscZcName << endl;

	os << "NUMDATES: " << numDates << endl;

	// Print the two set of index rates 
	for (idx=0; idx<numDates; idx++) {

		os << format(" %10s %10s %10s %10s %10s %15.5f `%s' ",
			DrlTDatePrint(NULL, mResetDates[idx]),
			DrlTDatePrint(NULL, mResetEffDates[idx]),
			DrlTDatePrint(NULL, mAccStartDates[idx]),
			DrlTDatePrint(NULL, mAccEndDates[idx]),
			DrlTDatePrint(NULL, mPayDates[idx]),
			mNotionals[idx],
			mFormulas[idx].c_str()); 
		os << endl;
		os << "Idx 1: " << (*mRateResets[idx]);
		os << endl;
		os << "Idx 2: " << (*mRateResets2[idx]);
		os << endl;
	}

	// Print dependencies
	this->KVPAtom::Put(os);


	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPFloatLeg2Idx::Get: failed for %s.\n", GetName());
    }
}



//---------------------------------------------------------------
// Only allowed two rate indices (with spread schedule),
// So only the two dependent rate indices in the first reset is printed.
// Other rate indices are assumed only differ by constant spreads.
//
ostream&
KVPFloatLeg2Idx::YacctionWriteRecursive(ostream& os)
{
static	char	routine[] = "KVPFloatLeg2Idx::YacctionWriteRecursive";
	int	idx;

    try {

	// Only two indices are allowed 
	if ( NumDep() > 2 )
		throw KFailure("%s: different floating rates at each reset "
			       "are not allowed in the wrapper.\n",
				routine);

	for (idx=0; idx<NumDep(); idx++)
        	Dep(idx)->YacctionWriteRecursive(os);

        YacctionWrite(os);

        return(os);

    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, GetName());
    }
}



//---------------------------------------------------------------
// Only allowed two rate indices (with spread schedule), 
// and one payment formula
//

ostream&
KVPFloatLeg2Idx::YacctionWrite(ostream& os, int indent)
{
	int	idx, numDates;
	bool	isOneRate = true;

    try {

	if (GetWriteFlag())
	{

	    numDates = mResetDates.size();

        CheckValid();

	    // Write the floating rate first
	    // Test if all the rates are the same for each reset
	    //
	    for (idx=1; idx<numDates; idx++) {
		if (!(mRateResets[0]->Rate()  == mRateResets[idx]->Rate()) &&
		    !(mRateResets2[0]->Rate() == mRateResets2[idx]->Rate()))
		isOneRate = false;
	    }

	    if (!isOneRate)
		throw KFailure("KVPFloatLeg2Idx::YacctionWrite: different "
			"floating rates are not allowed in the wrapper.\n");

	    os << GetName() << "=FLOATLEG({" << endl;
	    for (idx=0; idx<numDates; idx++) {

		os << format(" %10s %10s %10s %10s %10s",
			GtoFormatDate(mResetDates[idx]),
			GtoFormatDate(mResetEffDates[idx]),
			GtoFormatDate(mAccStartDates[idx]),
			GtoFormatDate(mAccEndDates[idx]),
			GtoFormatDate(mPayDates[idx]));
		
		os << " \"" << (mFormulas[idx].empty() ? 
			String("x0") : mFormulas[idx]) << "\" ";

		os << formatDouble("%15.5f", mNotionals[idx]);

		os << endl;
	    }


        if (!mDayCc.isNegative())
        {
            os << "}, "
               << DrlTDayCountPrint(NULL, (long)mDayCc) << ", "
               << "STUB_" << mStubConv << ", " << endl;
        }
        else  // risky dcc
        {
            os << "}, "
               << DrlTDayCountPrint(NULL, abs((int)mDayCc)) << "D, "
               << "STUB_" << mStubConv << ", " << endl;
        }
	
            // In case of same index
            if (NumDep() == 1)
            {
		if (Dep(0)->IsType("KRate"))
			os << "\t" << Dep(0)->GetName() << ", " 
                           << Dep(0)->GetName() << ", " << endl;
            }
            else   // 2 different indices
	        for (idx=0; idx<NumDep(); idx++)
	        {
		    if (Dep(idx)->IsType("KRate"))
			os << "\t" << Dep(idx)->GetName() << ", " << endl;
                }
	
	    /* replaced this with array of formulas */
	    // os << "\t\"" << (mFormulas[0].empty() ? 
	    //	String("nil") : mFormulas[0]) << "\","
	    //   << endl;

	    os << "\t\"" << mDiscZcName << "\");"
	       << endl << endl;

	    WriteDone();
	}
	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPFloatLeg2Idx::YacctionWrite: failed for %s.\n", GetName());
    }
}

