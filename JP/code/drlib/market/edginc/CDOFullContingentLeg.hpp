//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//
//----------------------------------------------------------------------------

#ifndef QR_CDOFULLCONTINGENTLEG_HPP
#define QR_CDOFULLCONTINGENTLEG_HPP
#include "edginc/CreditContingentLegBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** 'Full' contingent leg ie an implementation of ICDOContingentLeg which
    allows payAsYouGo and delayDays to vary per period */
class MARKET_DLL CDOFullContingentLeg: public CreditContingentLegBase{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    virtual ~CDOFullContingentLeg();

    /** Called immediately after object constructed. Overridden to validate
        length of delayDays and payAsYouGo arrays */
    virtual void validatePop2Object();

protected:
    /** Returns the number of delay days for payment for the specified
        period given by the supplied integer when 'paying as you
        go' (otherwise not defined). */
    virtual int numDelayDays(int periodIdx) const;

    /** Returns whether for the specified period we are 'paying as you
        go' */
    virtual bool payingAsYouGo(int periodIdx) const;

private:

    CDOFullContingentLeg();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    //-------
    // FIELDS
    //-------
    IntArray  delayDays;  // Number of delay days for payment when PorG='G'
    BoolArray payAsYouGo; // pay as you go (TRUE) or at payment date (FALSE)

};

DRLIB_END_NAMESPACE
#endif
