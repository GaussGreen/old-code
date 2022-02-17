//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InstrumentSettlement.hpp
//
//   Description : Defines how instruments can be settled
//
//   Author      : Andrew J Swain
//
//   Date        : 27 February 2001
//
//
//----------------------------------------------------------------------------


#ifndef INSTRUMENTSETTLEMENT_HPP
#define INSTRUMENTSETTLEMENT_HPP

#include "edginc/DateTime.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Asset.hpp"
#include "edginc/GetMarket.hpp"

DRLIB_BEGIN_NAMESPACE

/** There are many ways of settling an instrument - as cash, physically
    or as a margin option. This defines a set of interfaces that any
    implementation of instrument settlement need to obey.
*/

class MARKET_DLL InstrumentSettlement: public CObject,
                            virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;
    friend class InstrumentSettlementHelper;

    virtual ~InstrumentSettlement();

    /** given a trade date, when does this settle ? */
    virtual DateTime settles(const DateTime& date,
                             const CAsset*   asset) const = 0; // asset optional for N-factors
   
    /** is this a physical settlement ? */
    virtual bool isPhysical() const = 0;

    /** is this a margin option ? */
    virtual bool isMargin() const = 0;

    /** Compute discount factor between value date and the
        maturity date for the given date
    */
    virtual double pv(const DateTime& date, 
                      const YieldCurve* yc,
                      const Asset*    asset) const; // asset optional for N-factors

    /** Compute discount factor between value date and the
        maturity date for the given date
    */
    virtual double pv(const DateTime&   lodate,
                      const DateTime&   hidate,
                      const YieldCurve* yc,
                      const Asset*      asset) const; // asset optional for N-factors

    /** Compute discount factor between given date and the
        settlement date for the given date
    */
    virtual double pvAdjust(const DateTime& date, 
                            const YieldCurve* yc,
                            const Asset*    asset) const; // asset optional for N-factors

    /** populate from market cache - default implementation provided */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
protected:
    InstrumentSettlement(CClassConstSP clazz);
};

typedef smartConstPtr<InstrumentSettlement> InstrumentSettlementConstSP;
typedef smartPtr<InstrumentSettlement> InstrumentSettlementSP;
#ifndef QLIB_INSTRUMENTSETTLEMENT_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<InstrumentSettlement>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<InstrumentSettlement>);
EXTERN_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<InstrumentSettlementSP>(InstrumentSettlementSP* t));
EXTERN_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<InstrumentSettlementSP>(InstrumentSettlementSP* t, IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<InstrumentSettlement>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<InstrumentSettlement>);
INSTANTIATE_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<InstrumentSettlementSP>(InstrumentSettlementSP* t));
INSTANTIATE_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<InstrumentSettlementSP>(InstrumentSettlementSP* t, IObjectSP o));
#endif

DRLIB_END_NAMESPACE

#endif
