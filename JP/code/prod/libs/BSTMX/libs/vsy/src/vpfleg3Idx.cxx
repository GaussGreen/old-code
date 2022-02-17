/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      Robert Mattingly, David Liu
 ************************************************************************/
#include "vpfleg3Idx.h"

#include "kutilios.h"		// ios utilities


extern	"C" {
#include "datelist.h"           // GtoNewDateList 
#include "drltime.h"		// DrlTDatePrint

};



//---------------------------------------------------------------


KVPFloatLeg3Idx::KVPFloatLeg3Idx(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effective dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	const KVector(String)& formula,		// (I) payment formula array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const SharedPointer<KRate> floatRate1,  // (I) 1st rate index
	const SharedPointer<KRate> floatRate2,  // (I) 2nd rate index
	const SharedPointer<KRate> floatRate3,  // (I) 3rd rate index
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
static	char	routine[] = "KVPFloatLeg3Idx::KVPFloatLeg3Idx";
	int	numDates, idx;

	// Second index

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


	// Third index

	SharedPointer<KVPAtom>  rateIndex3;

	try {
	  
	  numDates = resetDates.size();
	  
	  for (idx=0; idx < numDates; idx++) {
	    mRateResets3.insert(mRateResets3.end(),
				new KRateReset(
					       resetDates[idx],
					       resetEffDates[idx],
					       *floatRate3)
				);
	  }
	  
	  
	  // Add rate as dependecy
	  SharedPointerConvertTo(floatRate3, rateIndex3);
	  AddDep(rateIndex3);
	}
	catch (KFailure) {
	  throw KFailure("%s: failed for %s.\n", routine, name);
	}
}



//---------------------------------------------------------------
//
KVPFloatLeg3Idx::~KVPFloatLeg3Idx()
{
	for(KVector(KRateReset*)::iterator it2=mRateResets2.begin();
                        it2!=mRateResets2.end(); ++it2)
                delete (*it2);
 
        mRateResets2.clear();

	for(KVector(KRateReset*)::iterator it3=mRateResets3.begin();
                        it3!=mRateResets3.end(); ++it3)
                delete (*it3);
 
        mRateResets3.clear();

}



//---------------------------------------------------------------
// Only called after construction

void
KVPFloatLeg3Idx::CheckValid() const
{
static	char	routine[] = "KVPFloatLeg3Idx::CheckValid";
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
	ASSERT_OR_THROW(mRateResets3.size()   == numDates);


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
	throw KFailure("KVPFloatLeg3Idx::CheckValid: failed for %s.\n", GetName());
    }

}


//---------------------------------------------------------------

ostream&
KVPFloatLeg3Idx::Put(ostream& os, int indent) const
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

	// Print the three sets of index rates 
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
		os << "Idx 3: " << (*mRateResets3[idx]);
		os << endl;
	}

	// Print dependencies
	this->KVPAtom::Put(os);


	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPFloatLeg3Idx::Put: failed for %s.\n", GetName());
    }
}



//---------------------------------------------------------------
// Only allowed three rate indices (with spread schedule),
// So only the two dependent rate indices in the first reset is printed.
// Other rate indices are assumed only differ by constant spreads.
//
ostream&
KVPFloatLeg3Idx::YacctionWriteRecursive(ostream& os)
{
static	char	routine[] = "KVPFloatLeg3Idx::YacctionWriteRecursive";
	int	idx;

    try {

	// Only three indices are allowed 
      if ( NumDep() > 3 )
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
KVPFloatLeg3Idx::YacctionWrite(ostream& os, int indent)
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
		    !(mRateResets2[0]->Rate() == mRateResets2[idx]->Rate())&&
		    !(mRateResets3[0]->Rate() == mRateResets3[idx]->Rate()))
		isOneRate = false;
	    }

	    if (!isOneRate)
		throw KFailure("KVPFloatLeg3Idx::YacctionWrite: different "
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
	throw KFailure("KVPFloatLeg3Idx::YacctionWrite: failed for %s.\n", GetName());
    }
}

