//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditPathValuesIn.hpp
//
//   Description : Generates instrument mtm values for an array of paths
//
//   Author      : Jay Blumenstein
//
//   Date        : 22 Aug 2002
//
//
//   $Log: CreditPathValuesIn.cpp,v $
//   Revision 1.9  2003/05/23 21:51:16  qhou
//   correct daily MTM infinite loop. daily MTM start from value date
//
//   Revision 1.8  2002/09/26 22:17:22  nshen
//   createDailySchedule() -> createMTMDates(). taking mtm_delay into account
//
//   Revision 1.7  2002/09/13 22:03:34  jblumens
//   added gemBandDates
//
//   Revision 1.6  2002/09/09 22:58:06  nshen
//   added registration for exposureCcy.
//
//   Revision 1.5  2002/09/09 17:56:12  nshen
//   added more fields.
//
//   Revision 1.4  2002/09/04 20:21:23  jblumens
//   a/2 +1 -> a/2, and some to do's
//
//   Revision 1.3  2002/09/03 17:51:11  nshen
//   make num of paths even if not - for anti path pairs.
//
//   Revision 1.2  2002/08/28 01:19:34  nshen
//   added createDailySchedule().
//
//   Revision 1.1  2002/08/27 20:25:58  jblumens
//   changed name,.
//
//   Revision 1.2  2002/08/26 23:17:56  nshen
//   a few additions.
//
//   Revision 1.1  2002/08/22 21:22:50  jblumens
//   registered some fields
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CreditEngine.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////
//
//  class CreditPathValuesIn
//
//  contains: All Credit inputs
//  methods: compute grids
//////////////////////////////////////////////////////////////////////

CClassConstSP const CreditPathValuesIn::TYPE = CClass::registerClassLoadMethod(
							"CreditPathValuesIn", typeid(CreditPathValuesIn), load);

void CreditPathValuesIn::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CreditPathValuesIn, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCreditPathValuesIn);
	FIELD(nPaths, "Number of paths in simulation");
	FIELD(mtmDates, "Dates to calc mark-to-market values on");
    FIELD_MAKE_OPTIONAL(mtmDates);
	FIELD(mtmDaily, "true = daily mtm (mtmDates not used), false = use mtmDates only");
	FIELD(exposureDates, "Dates to calculate default exposure on");
	FIELD(minCreditPPY, "Minimum periods per year for default");
	FIELD(positions, "sizes (with long/short). values are multiplied by positions");
	FIELD(seedForMarket, "seeds for market indexes");
	FIELD_MAKE_OPTIONAL(seedForMarket);
	FIELD(seedArr, "one seed per imnt");
	FIELD(correlationArr, "one corr for each (1-asset) imnt");
	FIELD_MAKE_OPTIONAL(correlationArr);
	FIELD(instArr, "array of instruments");
	FIELD(marketSeeds, "market seeds if more then one market");
	FIELD_MAKE_OPTIONAL(marketSeeds);
	FIELD(marketIndexes, "market index for each imnt");
	FIELD_MAKE_OPTIONAL(marketIndexes);
	FIELD(exposureCcy, "ccy in which exposure is calculated, default to USD.");
	FIELD_MAKE_OPTIONAL(exposureCcy);
    FIELD(gemBandDates, " dates for GEMS21 banding");
    FIELD(mtmDatesInUse, "Dates to calc mark-to-market values on -- after massage");   
    FIELD_MAKE_TRANSIENT(mtmDatesInUse);
	//    Addin::registerConstructor("CREDIT_Interface",
	//                                Addin::RISK,
	//                                "create input for Credit Interface",
	//                                TYPE);
}

void CreditPathValuesIn::validatePop2Object()
{
    const string method = "CreditPathValuesIn::validatePop2Object";

    if (instArr->size() == 0)
        throw ModelException(method, "empty instrument list.");
    if (seedArr.size() != instArr->size())
        throw ModelException(method, "num of seeds must equal to num of instruments.");

    // fill in market seeds if empty
    if( ! marketSeeds.size() )
    {
        if( seedForMarket == 0 )
            throw ModelException(method, "market seed is not specified or 0");

        marketSeeds.push_back( seedForMarket );
    }
    else
    {
        int n = marketSeeds.size();
        for( int i = 0; i < n; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                if( marketSeeds[ i ] == marketSeeds[ j ] )
                    throw ModelException(method, "market seeds must be different");
            }
        }
    }

    // fill in market indexes if empty
    if( ! marketIndexes.size() )
    {
        int n = instArr->size();
        marketIndexes.resize( n );
        for( int i = 0; i < n; ++i )
            marketIndexes[ i ] = 0;
    }
    else
    {
        if( marketIndexes.size() != instArr->size() )
            throw ModelException(method, "num of market indexes must equal to num of instruments.");

        int m = marketSeeds.size();
        int n = marketIndexes.size();
        for( int i = 0; i < n; ++i )
        {
            if( marketIndexes[ i ] < 0 || marketIndexes[ i ] >= m )
                throw ModelException(method, "market index is out of range");
        }
    }

    if (exposureDates.size() < 2)
        throw ModelException(method, "must have at least 2 exposure dates");
    // to do: more validations

    if (nPaths < 2)
        throw ModelException(method, "num of paths must be at least 2.");

    // always use even number for paths for antithetic pairs
    nPaths =2*(int)(nPaths/2);

    // to do:
	
	// check that no MTM or Default dates are in the past.
}

// if daily MTM, we need to generate dates starting from valueDate 
// due to path dependence of margin calls. 
void CreditPathValuesIn:: createMTMDates(HolidayConstSP hol, int mtm_delay, const DateTime valueDate)
{
    mtmDatesInUse.clear();

    if (mtmDaily)
    {
        // only need to consider mtm date prior to last exposure date - mtm_delay
        DateTime endDate = exposureDates[exposureDates.size()-1];
        
        DateTime aDate = hol->addBusinessDays(endDate, -mtm_delay-1);
        while (aDate >= valueDate )
        {
            mtmDatesInUse.push_back(aDate);
            aDate = hol->addBusinessDays(aDate,-1);
        }
    }
    else
    {
        mtmDatesInUse = mtmDates;
    }
}

DRLIB_END_NAMESPACE
