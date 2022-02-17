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

#ifndef QR_CDOCONTINGENTLEG_HPP
#define QR_CDOCONTINGENTLEG_HPP
#include "edginc/CreditContingentLegBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of ICDOContingentLeg which has all periods either as 'pay
    as you go' or 'pay at the end' */
class MARKET_DLL CDOContingentLeg: public CreditContingentLegBase {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    virtual ~CDOContingentLeg();
    
    /** Explicit constructor */
    CDOContingentLeg(
        int nbDelayDay,
        const DateTimeArray& obsStartDate,
        const DateTimeArray& obsEndDate,
        const DateTimeArray& payDate,
        const DoubleArray& notional,
        bool payAsYouGo);

protected:
    /** Returns the number of delay days for payment for the specified
        period given by the supplied integer when 'paying as you
        go' (otherwise not defined). The integer is ignored in this 
        implementation */
    virtual int numDelayDays(int periodIdx) const;

    /** Returns whether for the specified period we are 'paying as you
        go' The integer is ignored in this implementation */
    virtual bool payingAsYouGo(int periodIdx) const;

private:
    CDOContingentLeg();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    //-------
    // FIELDS
    //-------
    int           nbDelayDay; // Number of delay days for payment when PorG='G'
    bool          payAsYouGo; // pay at you go (TRUE) or at payment date (FALSE)

};

DRLIB_END_NAMESPACE
#endif
