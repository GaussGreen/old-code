//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LastSensDate.cpp
//
//   Description : Interface for instruments that know when to stop
//                 tweaking pointwise greeks
//
//   Author      : Andrew J Swain
//
//   Date        : 7 June 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/IModel.hpp"

DRLIB_BEGIN_NAMESPACE

LastSensDate::LastSensDate() {}

LastSensDate::~LastSensDate() {}

CClassConstSP const LastSensDate::TYPE = CClass::registerInterfaceLoadMethod(
    "LastSensDate", typeid(LastSensDate), 0);

// Attempt to delegate to product for endDate.
// This method may be used by instruments that implement the LastSensDate
// interface, but which delegate to the product.
DateTime LastSensDate::productEndDate(const CInstrument *instrument,
                                      const Sensitivity *sensControl) const {
    LastProductSensDate *lpsd = dynamic_cast<LastProductSensDate*>(
                                    sensControl->getModel());
    if (lpsd)
        return lpsd->endDate(instrument, sensControl);

    throw ModelException("LastSensDate::productSensDate",
        "A product-provided endDate is not implemented for model: " +
        sensControl->getModel()->getClass()->getName());
}

LastProductSensDate::LastProductSensDate() {}

LastProductSensDate::~LastProductSensDate() {}

CClassConstSP const LastProductSensDate::TYPE = CClass::registerInterfaceLoadMethod(
    "LastProductSensDate", typeid(LastProductSensDate), 0);

DRLIB_END_NAMESPACE
