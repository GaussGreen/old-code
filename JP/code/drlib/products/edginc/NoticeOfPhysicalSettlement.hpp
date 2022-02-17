//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Filename    : NoticeOfPhysicalSettlement.hpp
//
//   Description : Class used to store the "Notice of Physical Settlement", 
//                 ie, the amounts and ids of the obligations that will be 
//                 delivered as part of a default settlement.
//                 So far obligation ids are ignored
//
//   Author      : Jose Hilera
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_NOTICEOFPHYSICALSETTLEMENT_HPP
#define QLIB_NOTICEOFPHYSICALSETTLEMENT_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class NoticeOfPhysicalSettlement : public CObject {
public:
    static CClassConstSP const TYPE;
    ~NoticeOfPhysicalSettlement();
    
    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Returns the total % notional expected to arrive during this 
     * name's settlement period */
    double getNoPSNotionalFraction() const;

private:
    // For reflection
    static void load (CClassSP& clazz);
    NoticeOfPhysicalSettlement();
    static IObject* defaultNoticeOfPhysicalSettlement();
    
    // FIELDS
    CDoubleArray nopsAmounts; // Notional fraction amounts to be delivered
};

DECLARE(NoticeOfPhysicalSettlement);

DRLIB_END_NAMESPACE

#endif
