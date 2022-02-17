//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : FlatExpLossPrior.cpp
//
//   Description : Flat expected loss prior
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/FlatExpLossPrior.hpp"

DRLIB_BEGIN_NAMESPACE

// DEFAULT LOSS SET TO 0
const double FlatExpLossPrior::DEFAULT_LOSS = 0.0;

/** private constructor */
FlatExpLossPrior::FlatExpLossPrior() : CObject(TYPE) , expLoss(DEFAULT_LOSS)
{}

/** public constructor */
FlatExpLossPrior::FlatExpLossPrior(double expLoss) : CObject(TYPE), 
expLoss(expLoss)
{
	validatePop2Object();
}

/** Destructor */
FlatExpLossPrior::~FlatExpLossPrior()
{}

/** validate */
void FlatExpLossPrior::validatePop2Object()
{
	static const string method = "FlatExpLossPrior::validatePop2Object";
	try
	{
		if(expLoss <0 && expLoss > 1)
		{
			throw ModelException("expLoss (" + Format::toString(expLoss) +") is out of range [0,1]");
		}
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}

/** Returns expected loss surface for set of strikes and dates */
ExpectedLossSurfaceSP FlatExpLossPrior::getELSurface(
		const CDOQuotes & marketQuotes,			// market quotes
		const ICDSParSpreads & indexSwapSpreads	// index swap spreads
		) const
{
	static const string method = "FlatExpLossPrior::getELSurface";
	try
	{
		int t,s,i;
		/** strikes - just include 0 and 1 */
		DoubleArraySP strikes = DoubleArraySP(new DoubleArray(2));
		(*strikes)[0] = 0.0;
		(*strikes)[1] = 1.0;

		// define dates
		// use value date + indexSwapDates
		DateTime valueDate = marketQuotes.getValueDate();
		DateTimeArray valDateArray(1);
        valDateArray[0] = valueDate;
		
		/** index quote dates */ 
		int numIndexQuotes = indexSwapSpreads.getParSpreadsExpiries()->size();
		DateTimeArray indexQuoteDates(numIndexQuotes);
		 
		for(i = 0;i<numIndexQuotes;i++)
		{
			indexQuoteDates[i] = (*indexSwapSpreads.getParSpreadsExpiries())[i]->toDate(valueDate);
		}



		/** dates */
		DateTimeArraySP dates = DateTimeArraySP(new DateTimeArray(DateTime::merge(valDateArray,indexQuoteDates)));

		/** losses */
		DoubleArrayArraySP losses = DoubleArrayArraySP(new DoubleArrayArray(strikes->size()));
		
		
		for(s = 0; s < strikes->size(); s++)
		{
			(*losses)[s] = DoubleArray(dates->size());
			for(t = 0; t< dates->size(); t++)
			{
				(*losses)[s][t] = expLoss*(*strikes)[s]; // this gives base loss as fraction of portfolio notional =1
			}
		}


		return ExpectedLossSurfaceSP(new ExpectedLossSurface(dates,strikes,losses, indexSwapSpreads.getRecovery())); 
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}



void FlatExpLossPrior::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FlatExpLossPrior, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ITrancheQuoteInterpolator);
    EMPTY_SHELL_METHOD(defaultFlatExpLossPrior);

	FIELD(expLoss, "Constant expected loss as fraction of tranche size [default = 0]");
	FIELD_MAKE_OPTIONAL(expLoss);

}

IObject* FlatExpLossPrior::defaultFlatExpLossPrior() {
    return new FlatExpLossPrior();
}

CClassConstSP const FlatExpLossPrior::TYPE = 
    CClass::registerClassLoadMethod("FlatExpLossPrior", 
                                    typeid(FlatExpLossPrior), 
                                    load);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool FlatExpLossPriorLoad() {
    return (FlatExpLossPrior::TYPE != 0);
}


DRLIB_END_NAMESPACE
