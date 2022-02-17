//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Filename    : ClosedFormBSImpliedSmile.cpp
//
//   Description : Closed form pricing model for european options with
//                 a BS smile.
//
//   Author      : Charles Morcom
//
//   Date        : January 23, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormBSImpliedSmile.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketDataFetcherCDS.hpp"
#include "edginc/CDSVolCubeBSImpliedSmile.hpp"

DRLIB_BEGIN_NAMESPACE

ClosedFormBSImpliedSmile::~ClosedFormBSImpliedSmile() {}

ClosedFormBSImpliedSmile::ClosedFormBSImpliedSmile() : Model(TYPE) {}

ClosedFormBSImpliedSmile::ClosedFormBSImpliedSmile(CClassConstSP clazz): Model(clazz) {}


IModel::WantsRiskMapping ClosedFormBSImpliedSmile::wantsRiskMapping() const {
    return IModel::riskMappingIrrelevant;
}

/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP ClosedFormBSImpliedSmile::createMDF() const {
	// create MarketDataFetcher so that get CDS vols
  return MarketDataFetcherSP(new MarketDataFetcherCDS(false, true,CDSVolCubeBSImpliedSmile::TYPE));
}

/** calculate single price and store result in CResult */
void ClosedFormBSImpliedSmile::Price(CInstrument*  instrument, 
                            CControl*     control, 
                            CResults*     results){
    static const string method = "ClosedFormBSImpliedSmile::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormBSImpliedSmile::IntoProduct");
    }
    IProduct*     product = 0;
    try{
        product = intoProd->createProduct(this);
        product->price(this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}

//------------------------------
// IHasForwardRatePricer methods
//------------------------------

/** Key method providing access to the pricer */
IForwardRatePricerSP ClosedFormBSImpliedSmile::getForwardRatePricer() const
{
    if (!forwardRateModel)
    {
        //not provided by the user, so supply the default
        return getDefaultForwardRatePricer();
    }
    else
    {
        //use the supplied pricer
        return forwardRateModel;
    }
}
class ClosedFormBSImpliedSmileHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedFormBSImpliedSmile, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultClosedFormBSImpliedSmile);
        
    }

    static IObject* defaultClosedFormBSImpliedSmile() {
        return new ClosedFormBSImpliedSmile();
    }
};

CClassConstSP const ClosedFormBSImpliedSmile::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormBSImpliedSmile", typeid(ClosedFormBSImpliedSmile), ClosedFormBSImpliedSmileHelper::load);

// for linker
bool   ClosedFormBSImpliedSmileLoad() {
    return (ClosedFormBSImpliedSmile::TYPE != 0);
   }

CClassConstSP const ClosedFormBSImpliedSmile::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormBSImpliedSmile::IIntoProduct",
                                    typeid(ClosedFormBSImpliedSmile::IIntoProduct), 0);

DRLIB_END_NAMESPACE
