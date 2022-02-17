//----------------------------------------------------------------------------
//
//   File        : ClosedFormForwardRatePricer.cpp
//
//   Description : A closed form model capable
//                 of generating and pricing deterministic forward rates
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

ClosedFormForwardRatePricer::ClosedFormForwardRatePricer()
: CObject(TYPE)
{}

void ClosedFormForwardRatePricer::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(ClosedFormForwardRatePricer, clazz);
    SUPERCLASS(CObject); // ?? Model
    IMPLEMENTS(IForwardRatePricer);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

IObject* ClosedFormForwardRatePricer::defaultConstructor()
{
    return new ClosedFormForwardRatePricer();
}

CClassConstSP const ClosedFormForwardRatePricer::TYPE =
    CClass::registerClassLoadMethod(
        "ClosedFormForwardRatePricer",
        typeid(ClosedFormForwardRatePricer),
        load);

CClassConstSP const ClosedFormForwardRatePricer::IIntoProduct::TYPE =
    CClass::registerInterfaceLoadMethod(
        "ClosedFormForwardRatePricer::IIntoProduct",
        typeid(ClosedFormForwardRatePricer::IIntoProduct),
        0);

static void purifyFix(ClosedFormForwardRatePricer::IProduct  *product,
                      ClosedFormForwardRatePricer            *closeForm,
                      double                                 *result)
{
    *result = product->fwd(closeForm);
}


void ClosedFormForwardRatePricer::Forward(CObject *instrument, double *result)
{
    static const string method = "ClosedFormForwardRatePricer::Forward";
    IIntoProduct* intoProd;
    if(!IIntoProduct::TYPE->isInstance(instrument) ||
       !(intoProd = dynamic_cast<IIntoProduct*>(instrument)))
    {
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormForwardRatePricer::IntoProduct");
    }
    IProduct* product = 0;
    try{
        product = intoProd->createProduct(this);
        purifyFix(product, this, result);
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}


bool ClosedFormForwardRatePricerLoad() {
    return (ClosedFormForwardRatePricer::TYPE != 0);
}

DRLIB_END_NAMESPACE
