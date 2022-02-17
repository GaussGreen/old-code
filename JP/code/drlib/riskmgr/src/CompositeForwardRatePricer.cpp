//----------------------------------------------------------------------------
//
//   File        : CompositeForwardRatePricer.cpp
//
//   Description : A container for models capable
//                 of generating and pricing deterministic forward rates
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CompositeForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

CompositeForwardRatePricer::CompositeForwardRatePricer()
: CObject(TYPE)
{}

void CompositeForwardRatePricer::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CompositeForwardRatePricer, clazz);
    SUPERCLASS(CObject); // ?? Model
    IMPLEMENTS(IForwardRatePricer);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(irModel, "Model to be used for interest rate forwards");
    FIELD(crModel, "Model to be used for credit forwards");

    FIELD_MAKE_OPTIONAL(irModel); //one or other may be
    FIELD_MAKE_OPTIONAL(crModel); //specified
}

IObject* CompositeForwardRatePricer::defaultConstructor()
{
    return new CompositeForwardRatePricer();
}

void CompositeForwardRatePricer::Forward(CObject *instrument, double *result)
{
    static const string method = "CompositeForwardRatePricer::Forward";
    throw ModelException(method, " has not been implemented");
}



CClassConstSP const CompositeForwardRatePricer::TYPE =
    CClass::registerClassLoadMethod(
        "CompositeForwardRatePricer",
        typeid(CompositeForwardRatePricer),
        load);
/*
CClassConstSP const CompositeForwardRatePricer::IIntoProduct::TYPE =
    CClass::registerInterfaceLoadMethod(
        "CompositeForwardRatePricer::IIntoProduct",
        typeid(CompositeForwardRatePricer::IIntoProduct),
        0);
*/

bool CompositeForwardRatePricerLoad() {
    return (CompositeForwardRatePricer::TYPE != 0);
}

DRLIB_END_NAMESPACE
