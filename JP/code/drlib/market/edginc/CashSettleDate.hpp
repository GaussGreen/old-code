//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashSettleDate.hpp
//
//   Description : Cash settlement on a fixed date
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef CASHSETTLEDATE_HPP
#define CASHSETTLEDATE_HPP

#include "edginc/InstrumentSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

/** Cash settlement on a fixed date */

class MARKET_DLL CashSettleDate : public InstrumentSettlement {
public:
    static CClassConstSP const TYPE;

    CashSettleDate(const DateTime& settleDate);
    virtual ~CashSettleDate();
    
    /** given a trade date, when does this settle ? */
    virtual DateTime settles(const DateTime& date,
                             const CAsset*   asset) const;
   
    /** is this a physical settlement ? */
    virtual bool isPhysical() const;

    /** is this a margin option ? */
    virtual bool isMargin() const;

private:
    CashSettleDate();  // for reflection
    friend class CashSettleDateHelper;
    DateTime settleDate;
};

DRLIB_END_NAMESPACE

#endif
