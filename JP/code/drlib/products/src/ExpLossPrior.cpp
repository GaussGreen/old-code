//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : ExpLossPrior.cpp
//
//   Description : Custom expected loss prior
//
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/ExpLossPrior.hpp"

DRLIB_BEGIN_NAMESPACE


/** private constructor */
ExpLossPrior::ExpLossPrior() : CObject(TYPE) 
{}


/** Destructor */
ExpLossPrior::~ExpLossPrior()
{}

/** validate */
void ExpLossPrior::validatePop2Object()
{
	static const string method = "ExpLossPrior::validatePop2Object";
	try
	{
		int N = expiries->size();
		int M  = strikes->size();

		if(baseLosses->size() != M)
		{
			throw ModelException("BaseLosses array size ("+Format::toString(baseLosses->size())
				+") is not same as number of strikes ("+Format::toString(M)+")");
		}
		for(int i = 0;  i< M;i++)
		{
			if((*baseLosses)[i].size() != N)
			{
				throw ModelException("Size of doubleArray's in BaseLosses array  ("+Format::toString((*baseLosses)[i].size())
					+") is not same as number of expiries ("+Format::toString(N)+")");
			}
			for( int j = 0 ;j<N ; j++)
			{
				if((*baseLosses)[j][i] <0 && (*baseLosses)[j][i] > 1)
				{
					throw ModelException("baseLoss (" + Format::toString((*baseLosses)[j][i]) +") is out of range [0,1]");
				}
			}
			
		}
		
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}

/** Returns expected loss surface for set of strikes and dates */
ExpectedLossSurfaceSP ExpLossPrior::getELSurface(
	const CDOQuotes & marketQuotes,			// market quotes
	const ICDSParSpreads & indexSwapSpreads	// index swap spreads
	) const
{
	static const string method = "ExpLossPrior::getELSurface";
	try
	{
		int t;
		
		// define dates
		DateTime valueDate = marketQuotes.getValueDate();
		
		int N = expiries->size();
		/** dates */
		DateTimeArraySP dates = DateTimeArraySP(new DateTimeArray(N));
		for(t=0;t<N;t++)
		{
			(*dates)[t] = (*expiries)[t]->toDate(valueDate);
		}



		return ExpectedLossSurfaceSP(new ExpectedLossSurface(dates,strikes,baseLosses, indexSwapSpreads.getRecovery())); 
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}



void ExpLossPrior::load (CClassSP& clazz) {
	clazz->setPublic(); // make visible to EAS/spreadsheet
	REGISTER(ExpLossPrior, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(ITrancheQuoteInterpolator);
	EMPTY_SHELL_METHOD(defaultExpLossPrior);

	FIELD(expiries, "Expiries for losses");
	FIELD(strikes, "Base strikes for losses");
	FIELD(baseLosses, "Expected loss 'matrix'. Rows = dates, Cols = strikes");
	

}

IObject* ExpLossPrior::defaultExpLossPrior() {
	return new ExpLossPrior();
}

CClassConstSP const ExpLossPrior::TYPE = 
CClass::registerClassLoadMethod("ExpLossPrior", 
								typeid(ExpLossPrior), 
								load);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
* script to force the linker to include this file */
bool ExpLossPriorLoad() {
	return (ExpLossPrior::TYPE != 0);
}


DRLIB_END_NAMESPACE
