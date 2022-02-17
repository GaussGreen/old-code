//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ClosedFormCDSBasket.cpp
//
//   Description : Closed form model for pricing a basket of CDSs (plus some 
//                 basis adjustment, if required)
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormCDSBasket.hpp"
#include "edginc/MarketDataFetcherCIS.hpp"

DRLIB_BEGIN_NAMESPACE


ClosedFormCDSBasket::ClosedFormCDSBasket(): 
    CModel(TYPE), 
    calculateIndexBasis(true),
    useIndexBasis(true)
{}


/** Create the MarketDataFetcher used for retrieving market data etc */
MarketDataFetcherSP ClosedFormCDSBasket::createMDF() const {
    return MarketDataFetcherSP(
        new MarketDataFetcherCIS(true, calculateIndexBasis, useIndexBasis));
}

/** calculate single price and store result in CResult */
void ClosedFormCDSBasket::Price(CInstrument*  instrument, 
                                CControl*     control, 
                                CResults*     results)
{
    static const string method = "ClosedFormCDSBasket::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument)))
    {
        throw ModelException(method, 
                             "Instrument of type " +
                             instrument->getClass()->getName() +
                             " does not support ClosedFormCDSBasket::IIntoProduct");
    }

    IProduct* product = 0;
    try {
        product = intoProd->createProduct(this);
        product->price(this, control, results);
        //purifyFix(product, this, control, results);
    } catch (exception& e) {
        if (product) {
            delete product;
        }
        throw ModelException(e, method);
    }
    if (product) {
        delete product;
    }
}

IModel::WantsRiskMapping ClosedFormCDSBasket::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** Let the model determine which is the latest date it shows 
 * sensitivity to */
DateTime ClosedFormCDSBasket::endDate(const Sensitivity* sensitivity,
                                      const CInstrument* inst,
                                      const DateTime&    instEndDate) const 
{
    // This model does not do much yet (it is all mostly done in the instrument)
    // so return whatever the current last date is
    return instEndDate;
}

//------------------------------
// IHasForwardRatePricer methods
//------------------------------

/** Key method providing access to the pricer */
IForwardRatePricerSP ClosedFormCDSBasket::getForwardRatePricer() const
{
    if (!forwardRateModel)
    {
        //not provided by the user, so supply the default
        return getDefaultForwardRatePricer();
    }
    else
    {
        //use the supplied model
        return forwardRateModel;
    }
}

/** Invoked when Class is 'loaded' */
void ClosedFormCDSBasket::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ClosedFormCDSBasket, clazz);
    SUPERCLASS(Model);
    EMPTY_SHELL_METHOD(defaultClosedFormCDSBasket);
    FIELD(calculateIndexBasis, "If false, use the index basis provided "
                                      "in the market data. Otherwise (default) "
                                      "compute the index basis");
    FIELD(useIndexBasis,       "If false, no index basis adjustment will "
                                      "be used, ie, the instrument will be "
                                      "priced exactly as a basket of CDS on the "
                                      "underlying names. Default: true");
    FIELD       (forwardRateModel,    "A model capable of pricing all fees");

    FIELD_MAKE_OPTIONAL(calculateIndexBasis);
    FIELD_MAKE_OPTIONAL(useIndexBasis);
    FIELD_MAKE_OPTIONAL(forwardRateModel);
}


IObject* ClosedFormCDSBasket::defaultClosedFormCDSBasket() {
    return new ClosedFormCDSBasket();
}


CClassConstSP const ClosedFormCDSBasket::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormCDSBasket", 
    typeid(ClosedFormCDSBasket), 
    load);

CClassConstSP const ClosedFormCDSBasket::IIntoProduct::TYPE =
    CClass::registerInterfaceLoadMethod(
        "ClosedFormCDSBasket::IIntoProduct",
        typeid(ClosedFormCDSBasket::IIntoProduct), 0);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool ClosedFormCDSBasketLoad() {
    return (ClosedFormCDSBasket::TYPE != 0);
}

DRLIB_END_NAMESPACE
