//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ESWCashFlow.hpp
//
//   Description   cash flow class for ESW
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ESW_CASHFLOW_HPP
#define EDG_ESW_CASHFLOW_HPP

#include "edginc/Class.hpp"
#include "edginc/Results.hpp"
#include "edginc/OutputRequestUtil.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ESWCashFlow : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultESWCashFlow(){
        return new ESWCashFlow();
    }

    ESWCashFlow();

	/** validate CacsFlow leg */
	virtual void Validate();

	/** price cash flow leg */
	double priceCF(const DateTime& valDate, const DateTime& endDate,
        const DateTime& callSettleDate, bool isCallable,
		const Control* control, CResults* results,
        OutputRequestUtil::KnownCashFlows& knownCFs, 
        DateTimeArray& payment_dates) const;

	string		            LegName;	
    double                  SpotFX;
    DateTimeArray           PayDates;
    DoubleArray             Amounts;
    YieldCurveWrapper       discount;
};
typedef smartPtr<ESWCashFlow> ESWCashFlowSP;

DRLIB_END_NAMESPACE
#endif
