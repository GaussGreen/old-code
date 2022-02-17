//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericNFactor.hpp
//
//   Description : Base class for N factor generic instruments
//                 If you want to book an equity based N factor instrument
//                 in Pyramid 'Generics' this is the mandatory starting point
//
//   Author      : Andrew J Swain
//
//   Date        : 24 October 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_GENERICNFACTOR_HPP
#define EDR_GENERICNFACTOR_HPP

#include "edginc/Instrument.hpp"
// #include "edginc/InstrumentSettlement.hpp"
// #include "edginc/MultiMarketFactors.hpp"
#include "edginc/Theta.hpp"

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MarketData_forward.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(InstrumentSettlement);
FORWARD_DECLARE(IMultiMarketFactors);


/**  Base class for N factor generic instruments.
     If you want to book an equity based N factor instrument
     in Pyramid 'Generics' this is the mandatory starting point */
class PRODUCTS_DLL GenericNFactor: public CInstrument,
                      virtual public Theta::Shift {
public:
    static CClassConstSP const TYPE;

    virtual ~GenericNFactor();

    /** Implementation of abstract method in Instrument */
    virtual DateTime getValueDate()const;

    /** Validate instrument having aquired market data */
    void Validate();

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    //// roll through time (setting historic values)
    virtual bool sensShift(Theta* theta);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;


protected:
    GenericNFactor(CClassConstSP clazz);

    GenericNFactor(CClassConstSP clazz,
                   const DateTime&          valueDate,
                   double                   notional,         
                   InstrumentSettlementSP   instSettle, 
                   IMultiMarketFactorsSP    assets, 
                   const YieldCurveWrapper& discount);

    /** Adds the FWD_AT_MAT if requested */
    void addRequests(Control*        control,
                     Results*        results,
                     const DateTime& maturity) const;
    
    // fields
    DateTime                valueDate;
    double                  notional;         // Size of deal 
    InstrumentSettlementSP  instSettle;       // instrument settlement details 
    IMultiMarketFactorsSP   assets;           // the underlyings
    YieldCurveWrapper       discount;         // ccy to discount payoff 

private:
    friend class GenericNFactorHelper;

    GenericNFactor(); // not implemented
    GenericNFactor(const GenericNFactor& rhs);
    GenericNFactor& operator=(const GenericNFactor& rhs);
};

typedef smartConstPtr<GenericNFactor> GenericNFactorConstSP;
typedef smartPtr<GenericNFactor> GenericNFactorSP;

DRLIB_END_NAMESPACE
#endif

