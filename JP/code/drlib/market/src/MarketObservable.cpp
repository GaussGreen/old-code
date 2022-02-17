//----------------------------------------------------------------------------
//
//   Filename    : MarketObservable.cpp
//
//   Description : A market observable 
//                 This could be a price based asset (equity, commodity, fx)
//                 or a curve/strip of rates for energy/rates
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//   Author      : Ian Stares
//
//   Date        : January 31 2006
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MarketObservable.hpp"

DRLIB_BEGIN_NAMESPACE

void IMarketObservable::load(CClassSP& clazz)
{
    REGISTER_INTERFACE(IMarketObservable, clazz);
    EXTENDS(IObject);
}

const string IMarketObservable::DEFAULT_SOURCE = "DEFAULT";

ObservationSourceSP IMarketObservable::getDefaultObsSource() {
    return ObservationSourceSP(new ObservationSource(DEFAULT_SOURCE));
}

CClassConstSP const IMarketObservable::TYPE = CClass::registerInterfaceLoadMethod(
    "IMarketObservable", typeid(IMarketObservable), load);

DRLIB_END_NAMESPACE
