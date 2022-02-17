//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : vwap.hpp
//
//   Description : volume weighted average price instrument
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 26, 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VWAP_HPP
#define EDR_VWAP_HPP
#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MonteCarlo.hpp"



DRLIB_BEGIN_NAMESPACE

/** Shell structure for spot/hybrid/ratio payoffs */
class PRODUCTS_DLL VWAP: public CInstrument, 
                public Theta::IShift,
                public ISensitiveStrikes,
                public LastSensDate,
                public IMCIntoProduct  {
public:
    static CClassConstSP const TYPE;

    /** instrument validation */
    virtual void Validate();

    /** input data validation */
    virtual void validatePop2Object();

    /** retrieve market data needed by Vanilla - just valueDate, asset and
        discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** Implementation of MonteCarlo::IntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;


    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);
  
    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    friend class VWAPHelper;
    friend class VWAP_MC;

    friend class FastQuoteEnv;

    VWAP();
    VWAP(const VWAP& rhs);
    VWAP& operator=(const VWAP& rhs);
    static void load(CClassSP& clazz);

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

protected:
    VWAP(CClassConstSP clazz);

    DateTime                valueDate;
    DateTime                startDate;

    ScheduleSP              vwapSchedule;

    CAssetWrapper           asset;
    string                  ccyTreatment;
    YieldCurveWrapper       discount;
    InstrumentSettlementSP  instSettle;

    DateTime                lastDate;
    int                     nDays;
    HolidayWrapper          hols;
    DateTimeArray           sampleDates;
    
    double                  maxShares;
    double                  sharesSoFar;

    double                  maxDollars;
    double                  dollarsSoFar;

    DoubleArraySP           minPD;
    DoubleArraySP           maxPD;
    DoubleArraySP           shareMinMax;

    DoubleArraySP           minProfitPerShare;

    DoubleArraySP           histSharePrice;

    bool                    hasThreeDayRule;
    double                  threeDayLevel;
    int                     nDaysSoFar;

    int                     currentDateOffset;
};

typedef smartPtr<VWAP> VWAPSP;

DRLIB_END_NAMESPACE
#endif
