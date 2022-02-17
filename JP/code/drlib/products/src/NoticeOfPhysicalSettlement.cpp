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

#include "edginc/config.hpp"
#include "edginc/NoticeOfPhysicalSettlement.hpp"
#include "edginc/Format.hpp"


DRLIB_BEGIN_NAMESPACE


NoticeOfPhysicalSettlement::NoticeOfPhysicalSettlement() : CObject(TYPE)
{}

NoticeOfPhysicalSettlement::~NoticeOfPhysicalSettlement()
{}


/** Called immediately after object constructed */
void NoticeOfPhysicalSettlement::validatePop2Object() {
    static const string method = 
        "NoticeOfPhysicalSettlement::validatePop2Object";

    const int numOfNopsAmounts = nopsAmounts.size();

    if (numOfNopsAmounts == 0) {
        // Note we check for zero length rather than for 
        // getNoPSNotionalFraction() == 0  because a notional value of zero 
        // is valid: it can be used to indicate that no obligations will be 
        // delivered (and therefore no contingent payment will be made)
        throw ModelException(method,
                             "The NoPS amounts array cannot have zero lenght: "
                             "If the NoPS has not been delivered, do not pass "
                             "in the NoPS parameter - and if the NoPS indicate "
                             "that no obligations will be delivered, create "
                             "one entry in the array, with value 0.");
    }
    for (int i=0; i < nopsAmounts.size(); ++i) {
        if (nopsAmounts[i] < 0) {
            throw ModelException(method,
                                 "nopsAmounts[" + 
                                 Format::toString(i) + "] = " + 
                                 Format::toString(nopsAmounts[i]) + 
                                 ". Negative amounts are not allowed.");
        }
    }

    if (getNoPSNotionalFraction() > 1.0) { 
        throw ModelException(method,
                             "The sum of the notional fractions in the NoPS "
                             "is greater than 100% (" +
                             Format::toString(100.0 * getNoPSNotionalFraction()) +
                             ").");
    }
}


/** Returns the % notional expected to arrive during this 
 * name's settlement period */
double NoticeOfPhysicalSettlement::getNoPSNotionalFraction() const {
    double notionalFraction = 0.0;
    for (int i=0; i < nopsAmounts.size(); ++i) {
        notionalFraction += nopsAmounts[i];
    }
    return notionalFraction;
}


/** Invoked when Class is 'loaded' */
void NoticeOfPhysicalSettlement::load(CClassSP& clazz) {
    clazz->setPublic();  // make visible to EAS/spreadsheet
    REGISTER(NoticeOfPhysicalSettlement, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultNoticeOfPhysicalSettlement);

    FIELD(nopsAmounts, "Notice of Physical Delivery - Percentage of "
                              "notional amounts that will be delivered");
}

IObject* NoticeOfPhysicalSettlement::defaultNoticeOfPhysicalSettlement() {
    return new NoticeOfPhysicalSettlement();
}

CClassConstSP const NoticeOfPhysicalSettlement::TYPE = 
    CClass::registerClassLoadMethod("NoticeOfPhysicalSettlement", 
                                    typeid(NoticeOfPhysicalSettlement), 
                                    load);

/** Included in ProductsLib to force the linker to include this file */
bool NoticeOfPhysicalSettlementLoad() {
    return (NoticeOfPhysicalSettlement::TYPE != 0);
}


DRLIB_END_NAMESPACE
