//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Filename    : MultiQQuasiSmile.cpp
//
//   Description : Numerical pricing model for european options with
//                 a multiQ smile.
//
//   Author      : Mehdi Chaabouni
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MultiQQuasiSmile.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketDataFetcherCDS.hpp"

DRLIB_BEGIN_NAMESPACE

MultiQQuasiSmile::~MultiQQuasiSmile() {}

MultiQQuasiSmile::MultiQQuasiSmile() : Model(TYPE) {}

MultiQQuasiSmile::MultiQQuasiSmile(CClassConstSP clazz): Model(clazz) {}

IModel::WantsRiskMapping MultiQQuasiSmile::wantsRiskMapping() const {
    return IModel::riskMappingIrrelevant;
}

/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP MultiQQuasiSmile::createMDF() const {
    return MarketDataFetcherSP(new MarketDataFetcherCDS(true));
}

/** calculate single price and store result in CResult */
void MultiQQuasiSmile::Price(CInstrument*  instrument, 
                            CControl*     control, 
                            CResults*     results){
    static const string method = "MultiQQuasiSmile::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support MultiQQuasiSmile::IntoProduct");
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

//---------------------------------
// IHasForwardRatePricer methods
//---------------------------------
IForwardRatePricerSP MultiQQuasiSmile::getForwardRatePricer() const
{
    if(!forwardRateModel)
    {
        return getDefaultForwardRatePricer();
    }
    else
    {
        return forwardRateModel;
    }
}


class MultiQQuasiSmileHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MultiQQuasiSmile, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultMultiQQuasiSmile);
        FIELD (forwardRateModel, "A model capable of pricing all feees");
        FIELD_MAKE_OPTIONAL(forwardRateModel);
        
    }

    static IObject* defaultMultiQQuasiSmile() {
        return new MultiQQuasiSmile();
    }
};

CClassConstSP const MultiQQuasiSmile::TYPE = CClass::registerClassLoadMethod(
    "MultiQQuasiSmile", typeid(MultiQQuasiSmile), MultiQQuasiSmileHelper::load);

// for linker
bool   MultiQQuasiSmileLoad() {
    return (MultiQQuasiSmile::TYPE != 0);
   }

CClassConstSP const MultiQQuasiSmile::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("MultiQQuasiSmile::IIntoProduct",
                                    typeid(MultiQQuasiSmile::IIntoProduct), 0);

DRLIB_END_NAMESPACE
