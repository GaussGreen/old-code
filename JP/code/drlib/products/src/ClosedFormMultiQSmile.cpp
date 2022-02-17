//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Filename    : ClosedFormMultiQSmile.cpp
//
//   Description : Closed form pricing model for european options with
//                 a multiQ smile.
//
//   Author      : Charles Morcom
//
//   Date        : January 23, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormMultiQSmile.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketDataFetcherCDS.hpp"
#include "edginc/CDSVolCubeMultiQSmile.hpp"

DRLIB_BEGIN_NAMESPACE

ClosedFormMultiQSmile::~ClosedFormMultiQSmile() {}

ClosedFormMultiQSmile::ClosedFormMultiQSmile() : Model(TYPE) {}

ClosedFormMultiQSmile::ClosedFormMultiQSmile(CClassConstSP clazz): Model(clazz) {}

IModel::WantsRiskMapping ClosedFormMultiQSmile::wantsRiskMapping() const {
    return IModel::riskMappingIrrelevant;
}

/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP ClosedFormMultiQSmile::createMDF() const {
	// create fetcher so that get CDS vols but not IR vols
  return MarketDataFetcherSP(new MarketDataFetcherCDS(false,true,CDSVolCubeMultiQSmile::TYPE));
}

/** calculate single price and store result in CResult */
void ClosedFormMultiQSmile::Price(CInstrument*  instrument, 
                            CControl*     control, 
                            CResults*     results){
    static const string method = "ClosedFormMultiQSmile::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormMultiQSmile::IntoProduct");
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
IForwardRatePricerSP ClosedFormMultiQSmile::getForwardRatePricer() const
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
class ClosedFormMultiQSmileHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedFormMultiQSmile, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultClosedFormMultiQSmile);
        
    }

    static IObject* defaultClosedFormMultiQSmile() {
        return new ClosedFormMultiQSmile();
    }
};

CClassConstSP const ClosedFormMultiQSmile::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormMultiQSmile", typeid(ClosedFormMultiQSmile), ClosedFormMultiQSmileHelper::load);

// for linker
bool   ClosedFormMultiQSmileLoad() {
    return (ClosedFormMultiQSmile::TYPE != 0);
   }

CClassConstSP const ClosedFormMultiQSmile::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormMultiQSmile::IIntoProduct",
                                    typeid(ClosedFormMultiQSmile::IIntoProduct), 0);

DRLIB_END_NAMESPACE
