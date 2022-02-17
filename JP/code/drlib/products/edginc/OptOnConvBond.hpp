//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : OptOnConvBond.hpp
//
//   Description : Option on a Convertible Bond
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : April 25, 2002
//
//
//----------------------------------------------------------------------------

#ifndef OPTONCONVBOND_HPP
#define OPTONCONVBOND_HPP
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
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/ConvBond.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/IRiskyPricer.hpp"
#include "edginc/DerivativeAsset.hpp"
#include "edginc/IntrinsicMTM.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/FloatingBond.hpp"
#include "edginc/AssetSwapStrike.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL FixedYieldParameters: public AssetSwapStrike {
public:
    static CClassConstSP const TYPE;
    friend class OptOnConvBond;
    friend class FixedYieldParametersHelper;
    friend class FixedYieldAddin;

private:
    double          couponRate;
    int             couponFrequency;
    double          issuePrice;
    double          faceValue;
    DateTime        startDate;
    DateTime        endDate;
    double          yieldOverride;
    int             yieldType;
    HolidayWrapper  hols;

    FixedYieldParameters();

    double getYield() const;

    /** calculate the asset swap strike on a given date */
    virtual double getStrike(const DateTime&    baseDate,
                             ConvBondConstSP    cvb,
                             YieldCurveConstSP  discountCurve) const;

    /** returns true if the asset swap is exercisable on the baseDate */
    virtual bool   isExercisable(const DateTime& baseDate) const;

    /** reduced interface version of getStrike */
    virtual double getStrike(const DateTime&    baseDate) const;



    virtual void validatePop2Object(); // initialize a few params
};

typedef smartPtr<FixedYieldParameters>       FixedYieldParametersSP;
typedef smartConstPtr<FixedYieldParameters>  FixedYieldParametersConstSP;

//////// End of FixedYieldParameters Class  //////////////


class PRODUCTS_DLL LockOutParameters: public CObject {
public:
    static CClassConstSP const TYPE;
    friend class OptOnConvBond;
    friend class LockOutParametersHelper;

private:
    string            lockoutType;
    double            lockoutNotional;
    DateTime          lockoutDate;
    double            lockoutRate;
    string            lockoutDCCString;

    LockOutParameters();

    void validate();
};

typedef smartPtr<LockOutParameters>       LockOutParametersSP;
typedef smartConstPtr<LockOutParameters>  LockOutParametersConstSP;

/** Option on a Convertible Bond */

class PRODUCTS_DLL OptOnConvBond:         public CInstrument, 
                             public FD1F::IIntoProduct,
                             public FDModel::IIntoProduct,
                             public LastSensDate,
                             public ISensitiveStrikes,
                             public Theta::Shift,
                             public IScaleOutputs,
                             public IIntrinsicMTM,
                     virtual public ObjectIteration::IOverride, 
                             public IRiskyPricer {
public:
    static CClassConstSP const TYPE;
    friend class OptOnConvBondHelper;
    friend class OptOnConvBondClosedForm;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate(); // real validation

    virtual void validatePop2Object(); // initialize a few params

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual bool sensShift(Theta* shift);

    void   getMarket(const IModel* model, const MarketData* market);

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual double getFloatLegPV(DateTime  baseDate) const;

    virtual void getStrike(DateTime  baseDate, const double& bondFloor, bool includeLockout, 
                           bool *isExercisable, double *strike) const;
    
    virtual double lockoutFee(DateTime  baseDate) const;

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;
    
    /** Implementation of CFDGridPass::IntoProduct interface */
    virtual FD1F::IProduct* createProduct(FD1F* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    void recordOutputRequests(Control* control, Results* results, const double& convertPrice, 
                              const double& convertBondFloor, double theoOCBValue) const;

    /** breakdown of fixed leg calculation - used in addin to increase transparency 
        in pricing spreadsheet */
    ObjectArraySP calculateFloatingCashFlows(const DateTime& baseDate) const;

    /** breakdown of floating leg calculation - used in addin to increase transparency 
        in pricing spreadsheet */
    ObjectArraySP calculateFixedCashFlows(const DateTime& baseDate) const;

    ObjectArraySP calculateFixedAccrued(const DateTime& baseDate) const;

    // IRiskyPricer methods
    void setRisky(bool flag = true) { cvb->riskyGrowth = flag; } 
    bool isRisky() const { return cvb->riskyGrowth; }
    DateTime getEffMaturityDate() const { return cvb->effEndDate; }
    void setEffMaturityDate(const DateTime & maturityDate) {cvb->effEndDate = maturityDate; }

    // IScaleOutputs method
    void scaleOutputs(CControlSP control, ResultsSP unscaledResults);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    // IIntrinsicMTM method
    void calculateIntrinsic(CControlSP control, ResultsSP results);

    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

    // for product to access instrument data
    friend class ConvBondFDProd; 
    friend class ConvBond1fProd;

private:
    OptOnConvBond();
    OptOnConvBond(const OptOnConvBond& rhs);
    OptOnConvBond& operator=(const OptOnConvBond& rhs);

    double intrinsicValue(ResultsSP results) const;
    
    // -- Inputs -- // 
    DateTime                valueDate;
    CAssetWrapper           cvbAsset;

    string                  ocbType;
    ScheduleSP              exerSched;
    double                  spread;
    DateTime                swapMaturity;
    string                  swapIntervalString;
    bool                    shortFrontStub;
    double                  swapNotional;
    string                  swapDCCString;
    double                  swapBackEndFee;
    double                  swapLastFixRate;
    DateTime                swapLastFixDate;
    double                  swapBreakOfFundsRate;
    LockOutParametersSP     lockOutParameters;
    AssetSwapStrikeSP       assetSwapStrike;
    // FixedYieldParametersSP  fixedYieldParameters;

    // optional
    SettlementSP            bondSettle; // bond settlement: override the convertible bond settlement 

    // parameters for Asian unwinds
    bool                                      unwindNewSwap;
    FloatingBond::FloatingBondPaymentArraySP  swapFixings;

    // Derived
    ConvBondSP              cvb;
    DayCountConventionSP    swapDCC;
    MaturityPeriodSP        swapInterval;
    DayCountConventionSP    lockoutDCC;
    mutable double          bondFloor;
    double                  strikeAtOptionMat;
};


typedef smartConstPtr<OptOnConvBond> OptOnConvBondConstSP;
typedef smartPtr<OptOnConvBond> OptOnConvBondSP;


DRLIB_END_NAMESPACE
#endif
