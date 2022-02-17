/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      D. Liu
 ************************************************************************/
#include "vpfleg.h"

#include "kutilios.h"		// ios utilities


extern	"C" {
#include "datelist.h"       // GtoNewDateList 
#include "drltime.h"		// DrlTDatePrint

};


//---------------------------------------------------------------


KVPFloatLeg::KVPFloatLeg(
	const char *name,		// (I) name
	TDate startDate,		// (I) This date is not included
	TDate matDate,			// (I) maturity date
	const KDateInterval &freq,	// (I) frequency as interval
	KDayCc dayCc,			// (I) 
	KStubConv stubConv,		// (I) 
	TBoolean stubAtEnd,		// (I) T=stub at end, F=at beg 
	const SharedPointer<KRate> floatRate,	// (I) 
	const char *formula,		// (I)
	const char *discZcName)		// (I) 
	: KVPInstr(name)
{
static	char	routine[] = "KVPFloatLeg::KVPFloatLeg";
	TDateList	*datesDL = (TDateList *)NULL;
	int		idx;
	
	SharedPointer<KVPAtom>	rateIndex;

    try {
	// Set up reset effective dates. The first reset index is the next 
	// to last payment date (since we're moving backwards.)
	datesDL = GtoNewDateList(
		startDate,
		matDate, 
  		&((TDateInterval) freq),
		stubAtEnd);
	ASSERT_OR_THROW(datesDL != NULL);

	//
	// Check formula passed
	//
	if ((formula == NULL) ||
	    (!strlen(formula) != 0) ||
	    (!strcmp(formula, "NIL")) ||
	    (!strcmp(formula, "nil")) ||
	    (!strcmp(formula, "Nil")) ){
		formula = NULL;
	}


	// Assume accrueStarts=resets, and accrueEnds = payDates.
	//
	for (idx=0; idx < datesDL->fNumItems-1; idx++) {

		mResetDates.insert(mResetDates.end(), 
        			datesDL->fArray[idx]-floatRate->SpotOffset());
		mResetEffDates.insert(mResetEffDates.end(), 
        			datesDL->fArray[idx]);
		mAccStartDates.insert(mAccStartDates.end(), 
        				datesDL->fArray[idx]);
		mAccEndDates.insert(mAccEndDates.end(), 
        				datesDL->fArray[idx+1]);
		mPayDates.insert(mPayDates.end(), 
        				datesDL->fArray[idx+1]);
		mNotionals.insert(mNotionals.end(), 
        				1e0);

		mRateResets.insert(mRateResets.end(),
				   new KRateReset(
						mResetDates[idx],
						mResetEffDates[idx],
						*floatRate)
				   );

		mFormulas.insert(mFormulas.end(), 
        				String(formula ? formula : ""));
	}

	// Other info
	mDayCc = dayCc;
	mStubConv = stubConv;
	mDiscZcName = String(discZcName);



	// Add rate as dependecy
	SharedPointerConvertTo(floatRate, rateIndex);
	AddDep(rateIndex);

	GtoFreeDateList(datesDL);

    }
    catch (KFailure) {
	    GtoFreeDateList(datesDL);
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}


//---------------------------------------------------------------


KVPFloatLeg::KVPFloatLeg(
	const char *name,			// (I) object name
	const KVector(TDate)& resetEffDates,	// (I) reset effective dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const SharedPointer<KRate> floatRate,		// (I) paid rate
	const char *formula,			// (I) paymt formula (or NULL)
	const char *discZcName)			// (I) discount curve name
	: KVPInstr(name)
{
static	char	routine[] = "KVPFloatLeg::KVPFloatLeg";
	int		numDates, idx;

	SharedPointer<KVPAtom>  rateIndex;

    try {
	//
	// Test input data
	//
	if (resetEffDates.size() != accStartDates.size())
	    throw KFailure("%s: # resetEffDates (%d) != # accStartDates (%d).\n",
		routine, resetEffDates.size(), accStartDates.size());

	if (resetEffDates.size() != accEndDates.size())
	    throw KFailure("%s: # resetEffDates (%d) != # accEndDates (%d).\n",
		routine, resetEffDates.size(), accEndDates.size());

	if (resetEffDates.size() != payDates.size())
	    throw KFailure("%s: # resetEffDates (%d) != # payDates (%d).\n",
		routine, resetEffDates.size(), payDates.size());

	if (resetEffDates.size() != notionals.size())
	    throw KFailure("%s: # resetEffDates (%d) != # notionals (%d).\n",
		routine, resetEffDates.size(), notionals.size());

	numDates = resetEffDates.size();


	for (idx=0; idx < numDates; idx++) {
		if (accStartDates[idx] > accEndDates[idx])
			throw KFailure("%s: accStartDate %d (%s) >"
				" accEndDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accStartDates[idx]),
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]));
		if (accEndDates[idx] > payDates[idx])
			throw KFailure("%s: accEndDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

	}

	//
	// Check formula passed
	//
	if ((formula == NULL) ||
	    (!strlen(formula) != 0) ||
	    (!strcmp(formula, "NIL")) ||
	    (!strcmp(formula, "nil")) ||
	    (!strcmp(formula, "Nil"))) {
		formula = NULL;
	}

	//
	// Copy data
	//
	for (idx=0; idx < numDates; idx++) {
	    mResetDates.insert(   mResetDates.end(),
				  resetEffDates[idx]-floatRate->SpotOffset());
	    
	    if (mResetDates[idx] > payDates[idx])
			throw KFailure("%s: resetDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, mResetDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

	    mResetEffDates.insert(mResetEffDates.end(), resetEffDates[idx]);   
	    mAccStartDates.insert(mAccStartDates.end(), accStartDates[idx]);
	    mAccEndDates.insert(  mAccEndDates.end(),   accEndDates[idx]);
	    mPayDates.insert(     mPayDates.end(),      payDates[idx]);
	    mNotionals.insert(    mNotionals.end(),     notionals[idx]);
	    mFormulas.insert(     mFormulas.end(), 
        				String(formula ? formula : ""));

	    mRateResets.insert(mRateResets.end(),
			       new KRateReset(
					mResetDates[idx],
					mResetEffDates[idx],
					*floatRate)
			      );
	}

	// Other info
	mDayCc = dayCc;
	mStubConv = stubConv;
	mDiscZcName = String(discZcName);

	// Add rate as dependecy
	SharedPointerConvertTo(floatRate, rateIndex);
	AddDep(rateIndex);
    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}





//---------------------------------------------------------------


KVPFloatLeg::KVPFloatLeg(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effecitve dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const SharedPointer<KRate> floatRate,		// (I) paid rate
	const char *formula,			// (I) paymt formula (or NULL)
	const char *discZcName)			// (I) discount curve name
	: KVPInstr(name)
{
static	char	routine[] = "KVPFloatLeg::KVPFloatLeg";
	int		numDates, idx;

	SharedPointer<KVPAtom>  rateIndex;

    try {
	//
	// Test input data
	//
	if (resetDates.size() != resetEffDates.size())
	    throw KFailure("%s: # resetDates (%d) != # resetEffDates (%d).\n",
		routine, resetDates.size(), resetEffDates.size());

	if (resetDates.size() != accStartDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accStartDates (%d).\n",
		routine, resetDates.size(), accStartDates.size());

	if (resetDates.size() != accEndDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accEndDates (%d).\n",
		routine, resetDates.size(), accEndDates.size());

	if (resetDates.size() != payDates.size())
	    throw KFailure("%s: # resetDates (%d) != # payDates (%d).\n",
		routine, resetDates.size(), payDates.size());

	if (resetDates.size() != notionals.size())
	    throw KFailure("%s: # resetDates (%d) != # notionals (%d).\n",
		routine, resetDates.size(), notionals.size());

	numDates = resetDates.size();


	for (idx=0; idx < numDates; idx++) {
		if (accStartDates[idx] > accEndDates[idx])
			throw KFailure("%s: accStartDate %d (%s) >"
				" accEndDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accStartDates[idx]),
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]));

		if (resetDates[idx] > payDates[idx])
			throw KFailure("%s: resetDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, resetDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

		if (accEndDates[idx] > payDates[idx])
			throw KFailure("%s: accEndDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

	   	 mRateResets.insert(mRateResets.end(),
			            new KRateReset(
						resetDates[idx],
						resetEffDates[idx],
						*floatRate)
			      );
	}

	//
	// Check formula passed
	//
	if ((formula == NULL) ||
	    (!strlen(formula) != 0) ||
	    (!strcmp(formula, "NIL")) ||
	    (!strcmp(formula, "nil")) ||
	    (!strcmp(formula, "Nil"))) {
		formula = NULL;
	}

	//
	// Copy data
	//
	    mResetDates    = resetDates;
	    mResetEffDates = resetEffDates;   
	    mAccStartDates = accStartDates;
	    mAccEndDates   = accEndDates;
	    mPayDates      = payDates;
	    mNotionals     = notionals;
	    mFormulas.insert(mFormulas.end(), 
			     numDates, 
			     String(formula ? formula : ""));

	// Other info
	mDayCc = dayCc;
	mStubConv = stubConv;
	mDiscZcName = String(discZcName);

	// Add rate as dependecy
	SharedPointerConvertTo(floatRate, rateIndex);
	AddDep(rateIndex);
    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}



//---------------------------------------------------------------


KVPFloatLeg::KVPFloatLeg(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effecitve dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& spreads,		// (I) spread array
	const KVector(double)& notionals,	// (I) notional array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const SharedPointer<KRate> floatRate,		// (I) paid rate
	const char *formula,			// (I) paymt formula (or NULL)
	const char *discZcName)			// (I) discount curve name
	: KVPInstr(name)
{
static	char	routine[] = "KVPFloatLeg::KVPFloatLeg";
	int		numDates, idx;
	double		oldSprd, newSprd;

	SharedPointer<KVPAtom>  rateIndex;

	KRate		tmpRate(*floatRate);	// make a copy

    try {
	//
	// Test input data
	//
	if (resetDates.size() != resetEffDates.size())
	    throw KFailure("%s: # resetDates (%d) != # resetEffDates (%d).\n",
		routine, resetDates.size(), resetEffDates.size());

	if (resetDates.size() != accStartDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accStartDates (%d).\n",
		routine, resetDates.size(), accStartDates.size());

	if (resetDates.size() != accEndDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accEndDates (%d).\n",
		routine, resetDates.size(), accEndDates.size());

	if (resetDates.size() != payDates.size())
	    throw KFailure("%s: # resetDates (%d) != # payDates (%d).\n",
		routine, resetDates.size(), payDates.size());

	if (resetDates.size() != notionals.size())
	    throw KFailure("%s: # resetDates (%d) != # notionals (%d).\n",
		routine, resetDates.size(), notionals.size());

	if (resetDates.size() != spreads.size())
	    throw KFailure("%s: # resetDates (%d) != # spreads (%d).\n",
		routine, resetDates.size(), spreads.size());

	numDates = resetDates.size();


	// Original spread
	oldSprd = floatRate->Spread();

	for (idx=0; idx < numDates; idx++) {
		if (accStartDates[idx] > accEndDates[idx])
			throw KFailure("%s: accStartDate %d (%s) >"
				" accEndDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accStartDates[idx]),
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]));

		if (resetDates[idx] > payDates[idx])
			throw KFailure("%s: resetDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, resetDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

		if (accEndDates[idx] > payDates[idx])
			throw KFailure("%s: accEndDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

		// Add spread on top of current spread/coupon
		newSprd = oldSprd + spreads[idx];
		tmpRate.SetSpread(newSprd);

	   	 mRateResets.insert(mRateResets.end(),
			            new KRateReset(
						resetDates[idx],
						resetEffDates[idx],
						tmpRate)
			      );
	}

	//
	// Check formula passed
	//
	if ((formula == NULL) ||
	    (!strlen(formula) != 0) ||
	    (!strcmp(formula, "NIL")) ||
	    (!strcmp(formula, "nil")) ||
	    (!strcmp(formula, "Nil"))) {
		formula = NULL;
	}

	//
	// Copy data
	//
	    mResetDates    = resetDates;
	    mResetEffDates = resetEffDates;   
	    mAccStartDates = accStartDates;
	    mAccEndDates   = accEndDates;
	    mPayDates      = payDates;
	    mNotionals     = notionals;
	    mFormulas.insert(mFormulas.end(), 
			     numDates, 
			     String(formula ? formula : ""));

	// Other info
	mDayCc = dayCc;
	mStubConv = stubConv;
	mDiscZcName = String(discZcName);


	// Add rate as dependecy
	SharedPointerConvertTo(floatRate, rateIndex);
	AddDep(rateIndex);

    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}





//---------------------------------------------------------------


KVPFloatLeg::KVPFloatLeg(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effecitve dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	const KVector(String)& formulas,	// (I) pay formula array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const SharedPointer<KRate> floatRate,		// (I) paid rate
	const char *discZcName)			// (I) discount curve name
	: KVPInstr(name)
{
static	char	routine[] = "KVPFloatLeg::KVPFloatLeg";
	int		numDates, idx;

	SharedPointer<KVPAtom>  rateIndex;

    try {
	//
	// Test input data
	//
	if (resetDates.size() != resetEffDates.size())
	    throw KFailure("%s: # resetDates (%d) != # resetEffDates (%d).\n",
		routine, resetDates.size(), resetEffDates.size());

	if (resetDates.size() != accStartDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accStartDates (%d).\n",
		routine, resetDates.size(), accStartDates.size());

	if (resetDates.size() != accEndDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accEndDates (%d).\n",
		routine, resetDates.size(), accEndDates.size());

	if (resetDates.size() != payDates.size())
	    throw KFailure("%s: # resetDates (%d) != # payDates (%d).\n",
		routine, resetDates.size(), payDates.size());

	if (resetDates.size() != notionals.size())
	    throw KFailure("%s: # resetDates (%d) != # notionals (%d).\n",
		routine, resetDates.size(), notionals.size());

	if (resetDates.size() != formulas.size())
	    throw KFailure("%s: # resetDates (%d) != # formulas (%d).\n",
		routine, resetDates.size(), formulas.size());

	numDates = resetDates.size();

	mFormulas = formulas;
	for (idx=0; idx < numDates; idx++) {
		if (accStartDates[idx] > accEndDates[idx])
			throw KFailure("%s: accStartDate %d (%s) >"
				" accEndDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accStartDates[idx]),
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]));

		if (resetDates[idx] > payDates[idx])
			throw KFailure("%s: resetDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, resetDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

		if (accEndDates[idx] > payDates[idx])
			throw KFailure("%s: accEndDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

	   	 mRateResets.insert(mRateResets.end(),
			            new KRateReset(
						resetDates[idx],
						resetEffDates[idx],
						*floatRate)
			      );
		//
		// Check formula passed
		//
		if ((mFormulas[idx].empty())  ||
	    	    (mFormulas[idx] == "NIL") ||
	    	    (mFormulas[idx] == "nil") ||
	    	    (mFormulas[idx] == "Nil")) {
			mFormulas[idx] = "";
		}
	}


	//
	// Copy data
	//
	    mResetDates    = resetDates;
	    mResetEffDates = resetEffDates;   
	    mAccStartDates = accStartDates;
	    mAccEndDates   = accEndDates;
	    mPayDates      = payDates;
	    mNotionals     = notionals;

	// Other info
	mDayCc = dayCc;
	mStubConv = stubConv;
	mDiscZcName = String(discZcName);


	// Add rate as dependecy
	SharedPointerConvertTo(floatRate, rateIndex);
	AddDep(rateIndex);

    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}




//---------------------------------------------------------------

KVPFloatLeg::KVPFloatLeg(
	const char *name,			// (I) object name
	const KVector(TDate)& resetDates,	// (I) reset dates
	const KVector(TDate)& resetEffDates,	// (I) reset effecitve dates
	const KVector(TDate)& accStartDates,	// (I) acc start dates
	const KVector(TDate)& accEndDates,	// (I) acc end dates
	const KVector(TDate)& payDates,		// (I) pay dates array
	const KVector(double)& notionals,	// (I) notionals array
	const KVector(String)& formulas,	// (I) pay formula array
	const KVector(SharedPointer<KRate>) &floatRates, // (I) floating rate array
	KDayCc dayCc,				// (I) pay day count conv
	KStubConv stubConv,			// (I) stub conv
	const char *discZcName)			// (I) discount curve name
	: KVPInstr(name)
{
static	char	routine[] = "KVPFloatLeg::KVPFloatLeg";
	int		numDates, idx;

	SharedPointer<KVPAtom>  rateIndex;

    try {
	//
	// Test input data
	//
	if (resetDates.size() != resetEffDates.size())
	    throw KFailure("%s: # resetDates (%d) != # resetEffDates (%d).\n",
		routine, resetDates.size(), resetEffDates.size());

	if (resetDates.size() != accStartDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accStartDates (%d).\n",
		routine, resetDates.size(), accStartDates.size());

	if (resetDates.size() != accEndDates.size())
	    throw KFailure("%s: # resetDates (%d) != # accEndDates (%d).\n",
		routine, resetDates.size(), accEndDates.size());

	if (resetDates.size() != payDates.size())
	    throw KFailure("%s: # resetDates (%d) != # payDates (%d).\n",
		routine, resetDates.size(), payDates.size());

	if (resetDates.size() != notionals.size())
	    throw KFailure("%s: # resetDates (%d) != # notionals (%d).\n",
		routine, resetDates.size(), notionals.size());

	if (resetDates.size() != formulas.size())
	    throw KFailure("%s: # resetDates (%d) != # formulas (%d).\n",
		routine, resetDates.size(), formulas.size());

	if (resetDates.size() != floatRates.size())
	    throw KFailure("%s: # resetDates (%d) != # floatRates (%d).\n",
		routine, resetDates.size(), floatRates.size());


	numDates = resetDates.size();

	mFormulas = formulas;

	for (idx=0; idx < numDates; idx++) {
		if (accStartDates[idx] > accEndDates[idx])
			throw KFailure("%s: accStartDate %d (%s) >"
				" accEndDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accStartDates[idx]),
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]));

		if (accEndDates[idx] > payDates[idx])
			throw KFailure("%s: accEndDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, accEndDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));

		if (resetDates[idx] > payDates[idx])
			throw KFailure("%s: resetDate %d (%s) >"
				" payDate %d (%s).n", routine,
				idx+1, DrlTDatePrint(NULL, resetDates[idx]),
				idx+1, DrlTDatePrint(NULL, payDates[idx]));
		//
		// Check formula passed
		//
		if ((mFormulas[idx].empty())  ||
	    	    (mFormulas[idx] == "NIL") ||
	    	    (mFormulas[idx] == "nil") ||
	    	    (mFormulas[idx] == "Nil")) {
			mFormulas[idx] = "";
		}

		// Add rate as dependecy
		SharedPointerConvertTo(floatRates[idx], rateIndex);
		AddDep(rateIndex);

		// Insert rate array
		//
	   	mRateResets.insert(mRateResets.end(),
			            new KRateReset(
						resetDates[idx],
						resetEffDates[idx],
						*floatRates[idx]));
	}


	//
	// Copy data
	//
	mResetDates    = resetDates;
	mResetEffDates = resetEffDates;   
	mAccStartDates = accStartDates;
	mAccEndDates   = accEndDates;
	mPayDates      = payDates;
	mNotionals     = notionals;

	// Other info
	mDayCc = dayCc;
	mStubConv = stubConv;
	mDiscZcName = String(discZcName);

    }
    catch (KFailure) {
	throw KFailure("%s: failed for %s.\n", routine, name);
    }
}



//---------------------------------------------------------------
//
KVPFloatLeg::~KVPFloatLeg()
{
	for(KVector(KRateReset*)::iterator it=mRateResets.begin();
                        it!=mRateResets.end(); ++it)
                delete (*it);
 
	mRateResets.clear();

}



//---------------------------------------------------------------
// Only called after construction

void
KVPFloatLeg::CheckValid() const
{
static	char	routine[] = "KVPFloatLeg::CheckValid";
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
	throw KFailure("KVPFloatLeg::CheckValid: failed for %s.\n", GetName());
    }

}



//---------------------------------------------------------------

istream&
KVPFloatLeg::Get(istream& is, int drw)
{

    try {
	if (drw) {
	    int	idx, numDates;

	    mDayCc.Get(is, "KVPFloatLeg::Get: day count conv.");
	    mStubConv.Get(is, "KVPFloatLeg::Get: stub conv.");

	    numDates = getInt(is, "KVPFloatLeg::Get: number of dates.");

	    for (idx=0; idx<numDates; idx++) {
		mResetDates.insert(mResetDates.end(),
			getTDate(is, "KVPFloatLeg::Get: reset date."));
		mResetEffDates.insert(mResetEffDates.end(),
			getTDate(is, "KVPFloatLeg::Get: reset eff date."));
		mAccStartDates.insert(mAccStartDates.end(),
			getTDate(is, "KVPFloatLeg::Get: acc start date."));
		mAccEndDates.insert(mAccEndDates.end(),
			getTDate(is, "KVPFloatLeg::Get: acc end date."));
		mPayDates.insert(mPayDates.end(),
			getTDate(is, "KVPFloatLeg::Get: pay date."));
		mNotionals.insert(mNotionals.end(),
			getTDate(is, "KVPFloatLeg::Get: notional."));
	    }

	} else {
		throw KFailure("KVPFloatLeg::Get: format N/A.\n");
	}


	return(is);
    }
    catch (...) {
	throw KFailure("KVPFloatLeg::Get: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPFloatLeg::Put(ostream& os, int indent) const
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
	for (idx=0; idx<numDates; idx++) {

		os << format(" %10s %10s %10s %10s %10s %15.5f `%s' ",
			DrlTDatePrint(NULL, mResetDates[idx]),
			DrlTDatePrint(NULL, mResetEffDates[idx]),
			DrlTDatePrint(NULL, mAccStartDates[idx]),
			DrlTDatePrint(NULL, mAccEndDates[idx]),
			DrlTDatePrint(NULL, mPayDates[idx]),
			mNotionals[idx],
			mFormulas[idx].c_str());
		os << mRateResets[idx]->Rate();
		os << endl;
	}

	// Print dependencies
	this->KVPAtom::Put(os);


	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPFloatLeg::Get: failed for %s.\n", GetName());
    }
}



//---------------------------------------------------------------
// 
/*
void
KVPFloatLeg::SetWriteFlag(bool flag)
{
	int     idx, numDates;

	KVPAtom::SetWriteFlag(flag);
	
	// Set the rate flag
	numDates = mResetDates.size();
	
	for (idx=0; idx<numDates; idx++) {
		mResetDates[idx]->Rate().SetWriteFlag(flag);	
	}
		
}
*/



//---------------------------------------------------------------
// Only allowed one rate index (with spread schedule), 
// So only the first dependent rate index is printed.
// Other rate indices are assumed only differ by constant spreads.
// 
ostream&
KVPFloatLeg::YacctionWriteRecursive(ostream& os)
{
	Dep(0)->YacctionWriteRecursive(os);	
	
	YacctionWrite(os);
 
	return(os);
}




//---------------------------------------------------------------
// Only allowed one rate index (with spread schedule), 
// and one payment formula
//

ostream&
KVPFloatLeg::YacctionWrite(ostream& os, int indent) 
{
	int	idx, numDates;
	bool	isOneRate = true;
	double	spread;

    try {

	if (GetWriteFlag())
	{
 
	    numDates = mResetDates.size();

            CheckValid();

	    // Write the floating rate first
	    // Test if all the rates are the same for each reset
	    //
	    for (idx=1; idx<numDates; idx++) {
		if (!(mRateResets[0]->Rate() == mRateResets[idx]->Rate()) &&
		    mRateResets[idx]->Rate().IsFloating())
			isOneRate = false;
	    }

	    //
	    // YacctionWrite supports only single rate index in principle.
	    // Here is hack that would be allowed for step-up coupons
	    // or spreads, which are specified as spreads in the supyac
	    // schedule.
	    //
/*
	    if (!isOneRate)
		throw KFailure("KVPFloatLeg::YacctionWrite: different "
			"floating rates are not allowed in the wrapper.\n");
*/

	    os << GetName() << "=FLOATLEG({" << endl;
	    for (idx=0; idx<numDates; idx++) {

		// Spread depending on fixed/floating
		if (mRateResets[idx]->Rate().IsFloating())
			spread = mRateResets[idx]->Rate().Spread();
		else 	// fixed rate
			spread = mRateResets[idx]->Rate().Spread()
			       - mRateResets[0]->Rate().Spread();

		os << format(" %10s %10s %10s %10s %10s",
			GtoFormatDate(mResetDates[idx]),
			GtoFormatDate(mResetEffDates[idx]),
			GtoFormatDate(mAccStartDates[idx]),
			GtoFormatDate(mAccEndDates[idx]),
			GtoFormatDate(mPayDates[idx]));

		os << " \"" << (mFormulas[idx].empty() ? 
			String("x0") : mFormulas[idx]);

		if(!IS_ALMOST_ZERO(spread)){
		   os << "+" << formatDouble(" %10.5f", spread);
		}

		os << "\" ";

	        os << formatDouble("%12.2f", mNotionals[idx]);

		os << endl;
	    }

        if ((long)mDayCc > 0)
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
		break;
	    }
	
	    /* replaced this with array of formulas */
	    // os << "\t\"" << (mFormulas[0].empty() ? 
	    //	  String("nil") : mFormulas[0]) << "\","
	    //   << endl;

	    os << "\t\"" << mDiscZcName << "\");"
		<< endl << endl;

	    WriteDone();

	}

	return(os);
    }
    catch (KFailure) {
	throw KFailure("KVPFloatLeg::YacctionWrite: failed for %s.\n", GetName());
    }
}

