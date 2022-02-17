//----------------------------------------------------------------------------
//
//   File        : IForwardRatePricer.cpp
//
//   Description : A general interface for models that are capable
//                 of generating and pricing forward rates
//                 Also includes IHasForwardRatePricer interface
//                 that pricing models can implement to make their
//                 IForwardRatePricer accessible
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/Class.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

//---------------------------
// IForwardRatePricer methods
//---------------------------

const string IForwardRatePricer::INTEREST_RATE_FORWARD = "IR";
const string IForwardRatePricer::CREDIT_FORWARD = "CR";
const string IForwardRatePricer::CLOSED_FORM_FORWARD = "CF";

IForwardRatePricer::IForwardRatePricer(){}

IForwardRatePricer::~IForwardRatePricer(){}

void IForwardRatePricer::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IForwardRatePricer, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IForwardRatePricer::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IForwardRatePricer", typeid(IForwardRatePricer), load);

DECLARE(IForwardRatePricer)

//------------------------------
// IHasForwardRatePricer methods
//------------------------------

IHasForwardRatePricer::IHasForwardRatePricer(){}

IHasForwardRatePricer::~IHasForwardRatePricer(){}

IForwardRatePricerSP IHasForwardRatePricer::getDefaultForwardRatePricer() const
{
    ClosedFormForwardRatePricerSP pricer =
        ClosedFormForwardRatePricerSP(
            new ClosedFormForwardRatePricer());

    return pricer;
}

//does nothing
void IForwardRatePricer::getMarket(IModel *model, const MarketData *market)
{}


void IHasForwardRatePricer::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER_INTERFACE(IHasForwardRatePricer, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IHasForwardRatePricer::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IHasForwardRatePricer", typeid(IHasForwardRatePricer), load);

DECLARE(IHasForwardRatePricer)

//--------------------------------
// IWantsForwardRatePricer methods
//--------------------------------

IWantsForwardRatePricer::IWantsForwardRatePricer(){}

IWantsForwardRatePricer::~IWantsForwardRatePricer(){}

void IWantsForwardRatePricer::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER_INTERFACE(IWantsForwardRatePricer, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IWantsForwardRatePricer::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IWantsForwardRatePricer", typeid(IWantsForwardRatePricer), load);

DECLARE(IWantsForwardRatePricer)

DRLIB_END_NAMESPACE

