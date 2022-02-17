//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ClosedFormFA.cpp
//
//   Description : Closed form model for taking firm asset info as input
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 3, 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormFA.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/VolSurface.hpp"

DRLIB_BEGIN_NAMESPACE
ClosedFormFA::ClosedFormFA(): Model(TYPE), volType(VolSurface::TYPE->getName()), mdf(0) {}

static void purifyFix(ClosedFormFA::IProduct*     product,
                      ClosedFormFA*               closedFormFA,
                      CControl*                   control,
                      CResults*                   results){
    product->price(closedFormFA, control, results);
}

/** calculate single price and store result in CResult */
void ClosedFormFA::Price(CInstrument*  instrument, 
                         CControl*     control, 
                         CResults*     results){
    static const string method = "ClosedFormFA::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormFA::IntoProduct");
    }
    IProduct*     product = 0;
    try{
        product = intoProd->createProduct(this);
        //product->price(this, control, results);
        purifyFix(product, this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}


/** Override default createMDF in order to set the right MDF */
MarketDataFetcherSP ClosedFormFA::createMDF() const{
    return MarketDataFetcherSP(new MarketDataFetcherLN(volType));
}
    
IModel::WantsRiskMapping ClosedFormFA::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** Key method providing access to the pricer */
IForwardRatePricerSP ClosedFormFA::getForwardRatePricer() const
{
    return getDefaultForwardRatePricer();
}

class ClosedFormFAHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedFormFA, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultClosedFormFA);

        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
    }

    static IObject* defaultClosedFormFA(){
        return new ClosedFormFA();
    }
};

CClassConstSP const ClosedFormFA::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormFA", typeid(ClosedFormFA), ClosedFormFAHelper::load);
bool  ClosedFormFALoad() {
    return (ClosedFormFA::TYPE != 0);
   }


CClassConstSP const ClosedFormFA::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormFA::IIntoProduct",
                                    typeid(ClosedFormFA::IIntoProduct), 0);


DRLIB_END_NAMESPACE
