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
//----------------------------------------------------------------------------

#ifndef CREDIT_PATH_VALUE_GEN__H
#define CREDIT_PATH_VALUE_GEN__H

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Model.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////
//
//  class CreditPathValuesIn
//
//  contains: All inputs for generating values for each path at specified dates
//  methods: compute grids
//////////////////////////////////////////////////////////////////////
class CREDIT_DLL CreditPathValuesIn: public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
	
    static IObject* defaultCreditPathValuesIn(){
        return new CreditPathValuesIn();
    }

    virtual void validatePop2Object();

private:

    void createMTMDates(HolidayConstSP hol, int mtm_delay, const DateTime valueDate);

	CreditPathValuesIn() :  CObject(TYPE){seedForMarket = 0; exposureCcy="USD";}
    CreditPathValuesIn(const CreditPathValuesIn& rhs);
    CreditPathValuesIn& operator=(const CreditPathValuesIn& rhs);
	
    friend class CreditEngine;

    int				     nPaths;         /** Number of paths in simulation          */
    DateTimeArray	     mtmDates;       /** Dates to calc mark-to-market values on */
    bool                 mtmDaily;       /** true for daily mtm, false = read mtmDates */

    DateTimeArray	     exposureDates;   /** Dates to calculate default exposure on */
    int					 minCreditPPY;   /** Minimum periods per year for default   */

    int                  seedForMarket;     /** seed for market index			 */
	CIntArray			 seedArr;			/** one seed per imnt				 */
    DoubleArray			 correlationArr;    /** one corr for each (1-asset) imnt */
	CInstrumentArraySP	 instArr;			/** array of instruments */

    CIntArray            marketSeeds;     /** market seeds if more then one market */
    CIntArray            marketIndexes;   /** market index for each imnt */

    DoubleArray          positions;         /** sizes and long/short. path values are multiplied by positions */

    string               exposureCcy;       /**  exposure calculated in this ccy, fx fwd used to convert exposure values if needed */

    DateTimeArray        gemBandDates;       /** dates for GEMS21 banding */
    // below fields are not implemented at the moment
    // so not registered
    int					 curePeriod;        /** Cure period for unwinding position/hedge (in days).$unregistered */
    double  			 fireSaleDiscount;  /** fire sale discount $unregistered */
    bool                 longEquity;        /** JPM needs to sell equity on fire sale $unregistered */

    // transient
    DateTimeArray        mtmDatesInUse;
};

typedef smartPtr<CreditPathValuesIn> CreditPathValuesInSP;

DRLIB_END_NAMESPACE

#endif
