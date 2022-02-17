//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HandlePaymentEvents.cpp
//
//   Description : Interface that products can implement to indicate
//                 support for output requests for event handling
//                 e.g. PAYMENT_DATES & KNOWN_CASHFLOWS
//
//   Author      : Andrew J Swain
//
//   Date        : 22 August 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/HandlePaymentEvents.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IHandlePaymentEvents::TYPE = 
CClass::registerInterfaceLoadMethod("IHandlePaymentEvents", typeid(IHandlePaymentEvents), 0);

DRLIB_END_NAMESPACE
