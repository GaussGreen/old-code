//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaMoment.hpp
//
//   Description : Calculates moments (ATM vols, skew, convexity) of Vanilla Options
//
//   Author      : Regis Guichard
//
//   Date        : 27 March 03
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VANILLA_MOMENT_HPP
#define EDR_VANILLA_MOMENT_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/CompositeModel.hpp"
#include "edginc/Asset.hpp"
#include "edginc/InstrumentSettlement.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL VanillaMoment: public CInstrument,
                     public CompositeModel::IIntoProduct{
public:
    friend class VanillaMomentHelper;
    friend class VanillaMomentProduct;
    static CClassConstSP const TYPE;

    struct PRODUCTS_DLL RequestType{
        enum {
            ATM_VOL = 0,    // 0 + 1 strike  needed
            SKEW,           // 1 + 1 strikes needed
            CONVEXITY,      // 2 + 1 strikes needed
            NB_ENUMS
        };
    };
    
    virtual void validatePop2Object();

    /** instrument GetMarket to initiate market data selection.
        all it needs to do usually is to specify data tags for the market data
        this instrument depends on and then call data select in model.
        it allows instrument to implement specific get market actions */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    /** Called once before the initial pricing */
    virtual void Validate(){};

    /** Returns the value date (aka today) the instrument is currently
        pricing for */
    virtual DateTime getValueDate() const;

    virtual CompositeModel::IProduct* createProduct(CompositeModel* model) const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    static IObjectSP go(VanillaMoment* params);

    /** for reflection */
    VanillaMoment();

    DateTime                valueDate;
    CAssetWrapper           asset;
    string                  ccyTreatment;
    YieldCurveWrapper       discount;
    InstrumentSettlementSP  instSettle;
    InstrumentSettlementSP  premiumSettle;

    bool                    fwdStarting;
    DateTime                startDate;
    bool                    isCall;
    bool                    oneContract;
    double                  notional;
    double                  initialSpot;
    DateTimeArray           maturities;
    bool                    isAtmSpot;

    StringArray             requests;
};

DRLIB_END_NAMESPACE
#endif
