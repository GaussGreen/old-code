//----------------------------------------------------------------------------
//
//   File        : MQCMSPRICER.cpp
//
//   Description : A closed form model capable
//                 of generating and pricing deterministic forward rates
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MQCMSPricer.hpp"
#include "edginc/MQQuasiIRVolCMS.hpp"

DRLIB_BEGIN_NAMESPACE


MQCMSPricer::MQCMSPricer() : CObject(TYPE)
{
    isCashSettled = false;
    beta = 0.05;
    RichardsonOrder = 4;
}

void MQCMSPricer::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(MQCMSPricer, clazz);
    SUPERCLASS(CObject); // ?? Model
    IMPLEMENTS(IForwardRatePricer);

    FIELD(mqVol, "MQQuasiIRVolCMS, used for CMS adjustment");
    FIELD(isCashSettled, "default FALSE");
    FIELD(beta, "mean reversion, default 0.5");
    FIELD(RichardsonOrder, "integration order, default 4");

    FIELD_MAKE_OPTIONAL(isCashSettled);
    FIELD_MAKE_OPTIONAL(beta);
    FIELD_MAKE_OPTIONAL(RichardsonOrder);

    EMPTY_SHELL_METHOD(defaultConstructor);
}

IObject* MQCMSPricer::defaultConstructor()
{
    return new MQCMSPricer();
}

//get mqVol
MQQuasiIRVolCMSSP MQCMSPricer::getMQVol()
{
    return mqVol.getSP();
}
double MQCMSPricer::getBeta() const
{
    return beta;
}

bool MQCMSPricer::getIsCashSettled() const
{
    return isCashSettled;
}

double MQCMSPricer::getRichardsonOrder() const
{
    return RichardsonOrder;
}

static void purifyFix(MQCMSPricer::IProduct  *product,
                      MQCMSPricer            *closeForm,
                      double                 *result)
{
    *result = product->fwd(closeForm);
}

//forward
void MQCMSPricer::Forward(CObject* instrument, double *fwd)
{
    static const string method = "MQCMSPricer::Forward";
    IIntoProduct* intoProd;
    if(!IIntoProduct::TYPE->isInstance(instrument) ||
       !(intoProd = dynamic_cast<IIntoProduct*>(instrument)))
    {
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support MQCMSPricer::IntoProduct");
    }
    IProduct* product = 0;
    try{
        product = intoProd->createProduct(this);
        //product->price(this, control, results);
        purifyFix(product, this, fwd);
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}

//getMarket
void MQCMSPricer::getMarket(IModel *model, const MarketData *market)
{
    static const string method = "MQCMSPricer::getMarket";
    try {
        mqVol.getData(model, market);
    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}

CClassConstSP const MQCMSPricer::TYPE =
    CClass::registerClassLoadMethod("MQCMSPricer",
                                    typeid(MQCMSPricer),
                                    load);

CClassConstSP const MQCMSPricer::IIntoProduct::TYPE =
    CClass::registerInterfaceLoadMethod("MQCMSPricer::IIntoProduct",
                                        typeid(MQCMSPricer::IIntoProduct),
                                        0);

bool MQCMSPricerLoad() {
    return (MQCMSPricer::TYPE != 0);
}

DRLIB_END_NAMESPACE
