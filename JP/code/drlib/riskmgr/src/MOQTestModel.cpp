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

#include "edginc/config.hpp"

#include "edginc/MOQTestModel.hpp"

#include "edginc/MOQTestMarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE

MOQTestModel::MOQTestModel():
    CModel(TYPE)
{}

/** calculate single price and store result in CResult */
void MOQTestModel::Price(CInstrument*  instrument,
                         CControl*     control,
                         CResults*     results)
{
    throw ModelException("MOQTestModel::Price",
                         "This model does not price!");
}

/**
 * Whether to enable RiskMapping when computing sensitivities for
 * instruments priced using this model
 */
IModel::WantsRiskMapping MOQTestModel::wantsRiskMapping() const
{
    return riskMappingIrrelevant;
}

/** Creates an (MOQTest) MDF */
MarketDataFetcherSP MOQTestModel::createMDF() const
{
    return MarketDataFetcherSP(new MOQTestMarketDataFetcher(qualifierChoice));
}

CClassConstSP const MOQTestModel::TYPE =
    CClass::registerClassLoadMethod("MOQTestModel",
                                    typeid(MOQTestModel),
                                    load);

void MOQTestModel::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MOQTestModel, clazz);
    SUPERCLASS(CModel);
    EMPTY_SHELL_METHOD(MOQTestModel::defaultMOQTestModel);

    FIELD(qualifierChoice, "Choice of Market Data Qualifier string");
}

IObject* MOQTestModel::defaultMOQTestModel()
{
    return new MOQTestModel;
}

bool MOQTestModelLoad()
{
    return MOQTestModel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
