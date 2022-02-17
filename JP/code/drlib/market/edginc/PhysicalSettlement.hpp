//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PhysicalSettlement.hpp
//
//   Description : Physical style instrument settlement
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef PHYSICALSETTLEMENT_HPP
#define PHYSICALSETTLEMENT_HPP

#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Settlement.hpp"

DRLIB_BEGIN_NAMESPACE

/** Physical style instrument settlement */

class MARKET_DLL PhysicalSettlement : public InstrumentSettlement {
public:
    static CClassConstSP const TYPE;

    PhysicalSettlement();
    virtual ~PhysicalSettlement();
    
    /** given a trade date, when does this settle ? */
    virtual DateTime settles(const DateTime& date,
                             const CAsset*   asset) const;
   
    /** is this a physical settlement ? */
    virtual bool isPhysical() const;

    /** is this a margin option ? */
    virtual bool isMargin() const;

private:
    friend class PhysicalSettlementHelper;
};

DRLIB_END_NAMESPACE

#endif
