//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AtMaturity.hpp
//
//   Description : SA style 'at maturity' settlement
//
//   Author      : Andrew J Swain
//
//   Date        : 14 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef ATMATURITY_HPP
#define ATMATURITY_HPP

#include "edginc/InstrumentSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

/** SA style 'at maturity' settlement */

class MARKET_DLL AtMaturity : public InstrumentSettlement {
public:
    static CClassConstSP const TYPE;

    AtMaturity();
    virtual ~AtMaturity();
    
    /** given a trade date, when does this settle ? */
    virtual DateTime settles(const DateTime& date,
                             const CAsset*   asset) const;
   
    /** is this a physical settlement ? */
    virtual bool isPhysical() const;

    /** is this a margin option ? */
    virtual bool isMargin() const;

    /** Compute discount factor between value date and the
        settlement date for the given date
    */
    virtual double pv(const DateTime& date, 
                      const YieldCurve* yc,
                      const Asset*    asset) const;

    /** Compute discount factor between lodate and the
        maturity date for hidate */
    virtual double pv(const DateTime&   lodate,
                      const DateTime&   hidate,
                      const YieldCurve* yc,
                      const Asset*      asset) const;


    /** Compute discount factor between given date and the
        settlement date for the given date */
    virtual double pvAdjust(const DateTime& date, 
                            const YieldCurve* yc,
                            const Asset*    asset) const;

private:
    friend class AtMaturityHelper;
};

DRLIB_END_NAMESPACE

#endif
