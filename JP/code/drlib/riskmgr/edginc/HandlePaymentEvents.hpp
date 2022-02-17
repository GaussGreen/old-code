//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HandlePaymentEvents.hpp
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

#ifndef _HANDLEPAYMENTEVENTS_HPP
#define _HANDLEPAYMENTEVENTS_HPP

DRLIB_BEGIN_NAMESPACE

class Control;
class Results;

class RISKMGR_DLL IHandlePaymentEvents {
public:
    static CClassConstSP const TYPE;

    /** store events (requested in control) in results */
    virtual void recordEvents(Control* control,
                              Results* results) = 0;
};

// typedef smartConstPtr<IHandlePaymentEvents> IHandlePaymentEventsConstSP;
// typedef smartPtr<IHandlePaymentEvents> IHandlePaymentEventsSP;

DRLIB_END_NAMESPACE

#endif
