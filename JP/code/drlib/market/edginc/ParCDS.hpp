//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ParCDS.hpp
//
//   Description : Basic CDS implementation for pricing par instruments
//
//   Author      : Gordon Stephens
//
//   Date        : 19 April 2005
//
//
//----------------------------------------------------------------------------

#ifndef PARCDS_HPP
#define PARCDS_HPP

#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/CDSParSpreads.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ParCDS: public CInstrument, 
              public virtual ClosedForm::IIntoProduct,
              public LastSensDate
{
public:
    friend class ParCDSClosedForm;

    //--------
    // IObject
    //--------

    static CClassConstSP const TYPE;

    void validatePop2Object();

    //------------
    // CInstrument
    //------------

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    /** what's today ? */
    virtual DateTime getValueDate() const;

    virtual void Validate();

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    //-------------------------
    // ClosedForm::IIntoProduct
    //-------------------------

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;

    //-------------
    // LastSensDate
    //-------------

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    //-------
    // ParCDS
    //-------
    ParCDS(const string&  discountName,
           const string&  cdsParSpreadsName,
           DateTime       effectiveDate,
           ExpirySP       benchmark,
           double         spread,
           int            frequency);

private:

    //--------
    // IObject
    //--------

    static void load(CClassSP& clazz);

    static IObject* defaultParCDS();

    //-------
    // ParCDS
    //-------

    ParCDS();

    //closed form pricing
    void price(Control* control, CResults* results) const;

    //fields
    DateTime               valueDate;
    YieldCurveWrapper      discount;         /* Ccy to discount payoff                  */
    DateTime               effectiveDate;    /* When the CDS starts                     */
    ExpirySP               benchmark;        /* The maturity of the par instrument      */
    ICDSParSpreadsWrapper  cdsParSpreads;    /* The par curve                           */
    double                 spread;           /* the spread from the cdsParSpreads curve */
    int                    frequency;        /* the frequency of fee payments           */
};

typedef smartPtr<ParCDS> ParCDSSP;

DRLIB_END_NAMESPACE
#endif
