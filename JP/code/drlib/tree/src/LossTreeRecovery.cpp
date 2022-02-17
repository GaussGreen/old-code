//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : LossTreeRecovery.hpp
//
//   Description : Recovery model for spread loss tree
//					recovery can be a function of number of defaults that have incurred so far
//					Note that this has to be independent of the time of default
//
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/LossTreeRecovery.hpp"

DRLIB_BEGIN_NAMESPACE

/** private constructor */
LossTreeRecovery::LossTreeRecovery() : 
CObject(TYPE),
recoveries(),
totalLosses(),
totalRecoveries(),
numNames(-1)
{
}

void LossTreeRecovery::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(LossTreeRecovery, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultLossTreeRecovery);

	FIELD(recoveries, 
		"Array of recovery rate values in order of default. Length needs to be number of names in portfolio. "
		"[default = 0]"); // have default value for backwds compatibility
	FIELD_MAKE_OPTIONAL(recoveries);


	FIELD(numNames,"");
	FIELD_MAKE_TRANSIENT(numNames);

	FIELD(totalLosses,"");
	FIELD_MAKE_TRANSIENT(totalLosses);

	FIELD(totalRecoveries,"");
	FIELD_MAKE_TRANSIENT(totalRecoveries);
}



/** Destructor */
LossTreeRecovery::~LossTreeRecovery()
{}

IObject* LossTreeRecovery::defaultLossTreeRecovery() {
    return new LossTreeRecovery();
}

CClassConstSP const LossTreeRecovery::TYPE = 
    CClass::registerClassLoadMethod("LossTreeRecovery", 
                                    typeid(LossTreeRecovery), 
                                    load);


/** Called immediately after object constructed */
void LossTreeRecovery::validatePop2Object()
{
	static const string method = "LossTreeRecovery::validatePop2Object";
	
	// check if have recoveries input
	if(!!recoveries)
	{
		numNames = recoveries->size();

		int i;
		for(i = 0; i< numNames; i++)
		{
			if((*recoveries)[i] > 1 || (*recoveries)[i] < 0)
			{
				throw ModelException("recoveries[" + Format::toString(i) +"] = " + 
					Format::toString((*recoveries)[i]) +" is out of bounds [0,1]", method);
			}
		}

		initialise();
	}
	
		
}

/** set up recovery and loss arrays */
void LossTreeRecovery::initialise()
{
	if(!recoveries)
	{
		throw ModelException("recoveries are NULL", "LossTreeRecovery::initialise");
	}
	
	totalLosses = DoubleArraySP(new DoubleArray(numNames+1));
	totalRecoveries = DoubleArraySP(new DoubleArray(numNames+1));

	// scale all  total losses and recoveries by numNames
	(*totalLosses)[0] = 0;
	(*totalRecoveries)[0] = 0;
	
	int i;
	for(i = 1 ; i <= numNames ; i++)
	{
		(*totalLosses)[i] = (*totalLosses)[i-1] + (1 - (*recoveries)[i-1])/ ((double)numNames) ;
		(*totalRecoveries)[i] = (*totalRecoveries)[i-1] + (*recoveries)[i-1] / ((double)numNames);
	}
	

}

/** set the number of names in the portfolio */
void LossTreeRecovery::setNumNames(int numberOfNames)
{
	if( numberOfNames != numNames && numNames > 0)
	{
		throw ModelException("number of names in portfolio and recoveries array are not equal", 
			"LossTreeRecovery::setNumNames");
	}

	// initialies to 0 if no input given for backwds compatibility
	if(numNames < 1)
	{
		numNames = numberOfNames;
		recoveries = DoubleArraySP(new DoubleArray(numNames,0));

		initialise();
	}
		
}

double LossTreeRecovery::recovery(
	int numberOfDefaults // number of defaults incurred up to date 
	) const
{
	static const string method = "LossTreeRecovery::recovery";
	if (numNames <= 0 )
	{
		throw ModelException("NumNames < 0 or not set");
	}
	if(numberOfDefaults > numNames || numberOfDefaults < 0)
	{
		throw ModelException("numberOfDefaults out of bounds", method);
	}
	return (*recoveries)[numberOfDefaults];
}

 /** returns the total loss incurred up to date given the number of defaults */
double LossTreeRecovery::totalLoss(
	int numberOfDefaults // number of defaults incurred up to date 
	) const
{
	static const string method = "LossTreeRecovery::totalLoss";
	if (numNames <= 0 )
	{
		throw ModelException("NumNames < 0 or not set");
	}
	if(numberOfDefaults > numNames || numberOfDefaults < 0)
	{
		throw ModelException("numberOfDefaults out of bounds", method);
	}
	return (*totalLosses)[numberOfDefaults];
}

/** returns the total recovery incurred up to date given the number of defaults */
double LossTreeRecovery::totalRecovery(
								   int numberOfDefaults // number of defaults incurred up to date 
								   ) const
{
	static const string method = "LossTreeRecovery::totalRecovery";
	if (numNames <= 0 )
	{
		throw ModelException("NumNames < 0 or not set");
	}
	if(numberOfDefaults > numNames || numberOfDefaults < 0)
	{
		throw ModelException("numberOfDefaults out of bounds", method);
	}
	return (*totalRecoveries)[numberOfDefaults];
}

	

/** script to force the linker to include this file */
bool LossTreeRecoveryLoad() {
    return (LossTreeRecovery::TYPE != 0);
}


DRLIB_END_NAMESPACE
