//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids Derivatives Research
//
//   Filename    : IZeroCurveFactory.cpp
//
//   Description : Defines interface for zero curve factories.
//
//   Author      : Gordon C Stephens
//
//   Date        : 18 January 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZeroCurveFactory.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/RiskyCDSCurve.hpp"


DRLIB_BEGIN_NAMESPACE


/*
 * Reflection support.
 */

IZeroCurveFactory::IZeroCurveFactory(CClassConstSP clazz): CObject(clazz), 
                                                           adjustDates(false)
{
	//empty
}

void IZeroCurveFactory::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(IZeroCurveFactory, clazz);
    SUPERCLASS(CObject);
    
    FIELD(adjustDates, "adjust dates for backwards compatibility");
    
    FIELD_MAKE_TRANSIENT(adjustDates);
}


CClassConstSP const IZeroCurveFactory::TYPE = 
CClass::registerClassLoadMethod("IZeroCurveFactory", 
                                typeid(IZeroCurveFactory), 
                                load);


void IZeroCurveFactory::setAdjustDates(bool adjust)
{
    adjustDates = adjust;
}


// convert to risky curve
IYieldCurve* IZeroCurveFactory::makeRiskyCurve(
    const BootstrappedYieldCurve& riskless, 
    const CashFlowArray&          defaultRates, 
    double                        recovery) const
{
    static const string method = "IZeroCurveFactory::makeRisky";

    try
    {
        // create clean spread curve
        DoubleArray rates(defaultRates.size());
        ExpiryArray expiries(defaultRates.size());

        for (int i = 0; i < defaultRates.size(); ++i)
        {
            rates[i] = defaultRates[i].amount;
            expiries[i] = ExpirySP(new BenchmarkDate(defaultRates[i].date));
        }

        // create risky yield curve from clean spread curve
        CleanSpreadCurve* cs = new CleanSpreadCurve(riskless.getSpotDate(), &expiries, &rates);
        return new RiskyCDSCurve("risky", riskless, *cs, recovery);
    }
    catch(exception& e )
    {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
