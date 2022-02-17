//------------------------------------------------------------------------------
//
//   Group       : QR - Core Analytics 
//
//   Description : Dummy model for testing retrieval of market data using
//                 market data qualifiers
//
//   Author      : Andrew Greene 
//
//   Date        : 13 September 2006
//
//------------------------------------------------------------------------------

#ifndef MOQ_TEST_MODEL_HPP
#define MOQ_TEST_MODEL_HPP

#include "edginc/Model.hpp"

#include "edginc/MarketObjectQualifierString.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL MOQTestModel : public CModel
{
public:
    static CClassConstSP const TYPE;

    MOQTestModel();

    /** calculate single price and store result in CResult */
    void Price(CInstrument* instrument,
               CControl*    control,
               CResults*    results);

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     */
    WantsRiskMapping wantsRiskMapping() const;

protected:
    /** Creates an (MOQTest) MDF */
    MarketDataFetcherSP createMDF() const;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultMOQTestModel();

    string qualifierChoice;
};

DECLARE(MOQTestModel);

DRLIB_END_NAMESPACE

#endif // MOQ_TEST_MODEL_HPP
