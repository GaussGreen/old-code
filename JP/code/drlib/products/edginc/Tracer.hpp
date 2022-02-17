//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : Tracer.hpp
//
//   Description : Trigger Activated Convertible Security
//
//   Author      : Andr�Segger
//
//   Date        : 31 May 2002
//
//
//----------------------------------------------------------------------------

#ifndef TRACER_HPP
#define TRACER_HPP
#include "edginc/Instrument.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Bond.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/Asset.hpp"
#include "edginc/FD1F.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/ImpliedYTM.hpp"
#include "edginc/ImpliedYTP.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/ClosedForm.hpp"

DRLIB_BEGIN_NAMESPACE

/** Convertible Bond Instrument */

// TracerModel contains 2 models
// 1) to price a convertible bond
// 2) to price a double barrier
class PRODUCTS_DLL TracerModel: public CModel {
public:
    static CClassConstSP const TYPE;

    void validatePop2Object();

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(TracerModel* model,
                           Control*   control, 
                           Results*   results) const = 0;
        virtual ~IProduct(){};
    };
    
    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(const TracerModel* model) const = 0;
    };

    // inherited from CModel 
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    MarketObjectSP GetMarket(const MarketData*    market,
                             const string&        name,
                             const CClassConstSP& type) const;
    

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingDisallowed if either bond or barrier does,
     * else riskMappingAllowed if either bond or barrier does,
     * else riskMappingIrrelevant.
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;
    
    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz);

    // for TracerModel::IIntoProduct  
    static void loadIntoProduct(CClassSP& clazz);

    static IObject* defaultTracerModel();
    // constructor 
    TracerModel();

    // registered fields
    IModelSP    cvbModel;       // FD
    IModelSP    knockInModel;           // FD, Tree
    
protected:
    TracerModel(CClassConstSP clazz): CModel(clazz) {}

private:

};
typedef smartPtr<TracerModel> TracerModelSP;



class PRODUCTS_DLL Tracer: public CInstrument, 
              public LastSensDate,
              public Theta::Shift,
              public ISensitiveStrikes,
              public TracerModel::IIntoProduct {
public:
    static CClassConstSP const TYPE;
    friend class TracerHelper;
    friend class TracerProd;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate(); // real validation

    virtual void validatePop2Object(); // initialize a few params

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual bool sensShift(Theta* shift);

    /** create a tree1f payoff product */
    //virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;
    
    /** Implementation of TracerModel::IntoProduct interface */
    virtual TracerModel::IProduct* createProduct(const TracerModel* model) const;

    void recordOutputRequests(Control* control, Results* results, double fairValue) const;
    
    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    // for product to access instrument data
    friend class TracerFDProd;

private:
    Tracer();
    Tracer(const Tracer& rhs);
    Tracer& operator=(const Tracer& rhs);

    // instrument parameters
    YieldCurveWrapper         discount;
    CreditCurveWrapper        creditSpreads;
    CAssetWrapper             asset;
    DateTime                  valueDate;
    string                    ccyTreatment;

    // front bond parameters
    BondSP                    frontBond;
    double                    trigger;
    double                    conversionPremium;

    // back bond parameters - convertible
    MaturityPeriodSP          bondMaturity;
    string                    dayCountConvString;
    MaturityPeriodSP          callStartPeriod;
    double                    callTriggerLevel;
    double                    coupon;
    double                    redemption;
    double                    putLevel;
    MaturityPeriodSP          putStartPeriod;
    DateTimeArray             cbDates;
    int                       bondFrequency;
    HolidayWrapper            holidays;
    bool                      adjustRebateForAccrued;

    // back bond parameters - mandatory
    bool                      canConvIntoMandatory;
    double                    downsideProtection;
    double                    decsPremium;
    MaturityPeriodSP          decsMaturity;
    HolidayWrapper            decsHoliday;    
    double                    decsCoupon;    
    string                    decsDCC;
    int                       decsBondFrequency;
    bool                      decsDivPassThrough;

    bool                      backBondCapped;
    double                    backBondCap;

    bool                      knockedIn;
    DateTime                  knockInDate;

    //IModelSP                  cvbModel;
    //IModelSP                  knockInModel;

    // need a reference to the market to create the products
    CMarketDataSP             market;
};

typedef smartConstPtr<Tracer> TracerConstSP;
typedef smartPtr<Tracer>      TracerSP;

DRLIB_END_NAMESPACE
#endif
