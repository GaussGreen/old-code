//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CDSParSpreads.hpp
//
//   Description : Holds the current par spreads for CDSs
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : August 30, 2001
//
//
//----------------------------------------------------------------------------

#ifndef CDSPARSPREADS_HPP
#define CDSPARSPREADS_HPP

#include "edginc/MarketObject.hpp"

DRLIB_BEGIN_NAMESPACE
class CDSParSpreads;
#ifndef QLIB_CDSPARSPREADS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CDSParSpreads>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CDSParSpreads>);
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CDSParSpreads>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CDSParSpreads>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CDSParSpreads>);
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CDSParSpreads>);
#endif

DRLIB_END_NAMESPACE

#include "edginc/CDSParSpreadsBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** Holds the current par spreads for CDSs */
class MARKET_DLL CDSParSpreads: public CDSParSpreadsBase,
                     public virtual CreditDefaultSensBase::IShift
{
public:
    static CClassConstSP const TYPE;
    friend class CDSParSpreadsHelper;

    virtual ~CDSParSpreads();

    /** A constructor for convenience to be called by other QLib classs */
    CDSParSpreads(string                   name,
                  DateTime                 valueDate,
                  int                      spotOffset,
                  double                   parRecovery,
                  int                      parSwapFreq,
                  ExpiryArraySP            expiries,
                  CDoubleArraySP           spreads,
                  CDoubleArraySP           upfronts,
                  bool                     parAccrueFee,
                  string                   parDCC,       // day count convention for par CDS
                  string                   parBDC,       // bad day convention for par CDS
                  HolidayConstSP           parHols,      // holidays for par CDS
                  YieldCurveConstSP        discount,     // corresponding discount curve
                  DecretionCurveConstSP    prepay        // prepay curve
                  );
    
    
    /** Adjust dates */
    CDSParSpreadsSP adjustDates(string name, ExpiryArraySP adjEndDates) const;

    /** Get the par spreads curve from the market data cache */
    virtual void getMarket(const IModel* model, const MarketData *market);
    
    /** Returns name identifying CDS Par Curve for CREDIT_DEFAULT_SENS */
    virtual string sensName(CreditDefaultSensBase* sens) const;
    /** Shifts the object using given shift. */
    virtual bool sensShift(CreditDefaultSensBase* sens);

    /** Returns true if the name has defaulted */
    virtual bool defaulted() const;
    
    /** Returns the date that this name defaulted. If a default has not
        happened an 'empty' DateTime is returned */
    virtual const DateTime& getDefaultDate() const;

private:
    CDSParSpreads();

    static void acceptWrapperNameCollector(const CDSParSpreads* parSpreads,
                                           WrapperNameCollector* collector);

    /// non-static fields ////////
    bool                   hasDefaulted; // has a default happened
    DateTime               defaultDate;  // if so, when
};

typedef smartConstPtr<CDSParSpreads> CDSParSpreadsConstSP;
typedef smartPtr<CDSParSpreads>      CDSParSpreadsSP;
typedef MarketWrapper<CDSParSpreads> CDSParSpreadsWrapper;
DRLIB_END_NAMESPACE

#endif
