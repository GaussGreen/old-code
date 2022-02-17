//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DblBarrier.hpp
//
//   Description : double barrier instrument. 
//                 A KO barrier takes priority over KI barrier.
//                 If both are KI barriers, one KI counts as in.
//
//----------------------------------------------------------------------------

#ifndef EDG_DBL_BARRIER_HPP
#define EDG_DBL_BARRIER_HPP

#include "edginc/Generic1Factor.hpp"
#include "edginc/Asset.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/RiskyCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Theta.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/FD1F.hpp"
#include "edginc/FD1FDDE.hpp"
#include "edginc/ITaxableInst.hpp"
#include "edginc/LegalTerms.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

class CDblBarrier : public Generic1Factor, 
                    virtual public FD1F::IIntoProduct,
                    virtual public FD1FGeneric::IIntoProduct,
                    virtual public FDModel::IIntoProduct,
                    virtual public LastSensDate,
                    virtual public ITaxableInst::Basic,
                    virtual public LegalTerms::Shift        
{
public:
    static CClassConstSP const TYPE;

    // override base implementation if required
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    virtual void Validate();
 
    // below are copied from Vanilla
    virtual DateTime getValueDate() const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;
    
    bool sensShift(Theta* shift);


    /** create a fd payoff product */
    virtual FD1F::IProduct* createProduct(FD1F* model) const;

    /** create a generic fd payoff product */
    virtual FD1FGeneric::IProduct* createProduct(FD1FGeneric* model) const;

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel *     model);

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    CDblBarrier(const DateTime&           valueDate, 
                const DateTime&        startDate,
                bool                   fwdStarting,
                bool                   oneContract,
                double                 notional,
                double                 initialSpot,
                InstrumentSettlementSP instSettle,
                InstrumentSettlementSP premiumSettle,
                YieldCurveWrapper      discount,
                CAssetWrapper          asset,
                const string&          ccyTreatment,
                const string&          PayoffMode,
                ScheduleSP             exerciseSchedule,
                bool                   canExerciseEarly,
                ScheduleSP             upperBarrier,
                const string&          upperBarType,    
                ScheduleSP             upperRebate,
                ScheduleSP             lowerBarrier,
                const string&          lowerBarType,
                ScheduleSP             lowerRebate,
                bool                   intraDayMonitor,
                bool                   rebateAtMat,
                bool                   rebateNotScaled,
                const string&          barrierDependence,
                const string&          monitoringDependence);
 
    /** for ITaxableInst::Basic */
    const DateTime getFinalPaymentDate() const;

    /** for LegalTerms::Shift */
    virtual bool sensShift(LegalTerms* shift);

    /** add barrier dates to critical date list */
    static void addCritBarDates(const ScheduleSP bar, const DateTime& valueDate, const DateTime& matDate,
                            CTree1f* tree1f, DateTimeArray& critDates);

private:
    friend class CDblBarrierHelper;
    friend class CDblBarrier1fProd;
    friend class DblBarrierFDProd;

protected:
    static void load(CClassSP& clazz);
    CDblBarrier();
    CDblBarrier(CClassConstSP clazz);

    // this block is the same as Vanilla
    string              PayoffMode; // registered
    bool                isCall; // not registered
    bool                canExerciseEarly;
    ScheduleSP          exerciseSchedule;

    double              spotAtMaturity;

    DateTime            dateExercised;
    bool                isExercised;

    // barrier data
    ScheduleSP          UpperBarrier;
    ScheduleSP          LowerBarrier;

    string              UpperBarType;
    string              LowerBarType;

    bool                IntraDayMonitor; // false means once in a day

    ScheduleSP          UpperRebate;
    ScheduleSP          LowerRebate;
    bool                RebateAtMat;
    bool                RebateNotScaled;

    // above is RISK information i.e. what we price with
    // these are ECONOMIC information i.e. what we should report
    ScheduleSP          upperEcoBarrier;
    ScheduleSP          lowerEcoBarrier;
    
    // below four fields may be modified in priceDeadInstrument() const
    mutable bool        UpperBarBreached;
    mutable DateTime    UpperBarBreachDate;
    mutable bool        LowerBarBreached;
    mutable DateTime    LowerBarBreachDate;

    string              BarrierDependence;                 
    string              MonitoringDependence;                 

    int                 DEBUG_num_segments; //Number of segment. (for debug use).
    double              DEBUG_SmoothBarrierWidth; // negative value (default) is no smooth.  double value are used as SmoothWidth, e.g. 0.05 (5%) shift.

    // non registered.  Required as we need to know if after a tweak whether the barrier has been breached.
    mutable bool UpperBarBreachedDynamic; // $unregistered
    mutable bool LowerBarBreachedDynamic; // $unregistered
};

typedef smartPtr<CDblBarrier> CDblBarrierSP;

/////////////////////////////////////////////////////////
//           tree1f/fd product
/////////////////////////////////////////////////////////
/** DblBarrier product payoff for a tree */
// old interface and just keep it for now to allow FD1F and DDEd working
class CDblBarrier1fProd:  virtual public FD1F::IProduct,
virtual public FD1FGeneric::IProduct,  virtual public IDDEInitiator
{
public:
    friend class CDblBarrier; // in order to use isStaticSpread 
    
    typedef enum{CALL, PUT, BINARY, FORWARD} TPayoffMode;
    typedef enum{NA=-1, KI=0, KO=1} TBarrier;
    typedef enum{KI_KEEP_KO, KI_CANCEL_KO, ONCE_TOUCH, TWO_TOUCH} TBarDependence;
    typedef enum{BOTH, UPPER, LOWER} TMonDependence;

    CDblBarrier1fProd(const CDblBarrier* instr);
    
    virtual ~CDblBarrier1fProd(){}

    virtual CAssetConstSP GetAssetRef();

    virtual bool GetFwdStartLV();

    virtual DateTime GetFwdStartDateLV();

    virtual YieldCurveConstSP GetDiscCurveRef();

    virtual CVolRequestConstSP GetLNRequest();

    virtual void InitFD(CControl* control);

    /** initialise product specific data */
    virtual void InitProd();

    /** calculate at barriers for FD */
    virtual void preCalcFD(int step, int idx, int pStart, int pEnd);

    /** product payoff method at maturity */
    virtual void PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                             double * const * price);
    /** product payoff method at steps earlier than maturity */
    virtual void PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                 double * const * price);

    /** premium scaling */
    virtual double scalePremium(const double& fairValue);

    /** extra output requests */
    virtual void recordOutputRequests(Control* control, Results* results, double fairValue);

    /** extra output requests */
    virtual bool Positive() { return (PayMode != FORWARD); }

    virtual string getCcyTreatment();

    /** make price refinement - control variate */
    virtual double RefinePrice(double basePrice, double discFactor, bool) {
        // scale premium for forward starting.
        double price = discFactor*scalePremium(basePrice);
        return price;
    }
    
    virtual double getCoupon(int step, const double* s, int start, int end)
    { return 0.0; }
    
    virtual bool hasEquityLayer()
    { return false; }

    // for DDE
    DateTime maxMaturity() const;
    
    void sensitiveDates(  DateTimeArray    &dates) const;
    
    void sensitiveStrikes(  const DateTimeArray     dates,
                            DoubleArray             &strikes,   // same dimension as dates
                            bool                    &strikeIsPct) const;
       
protected:
 
    /** price dead instrument */
    bool    CheckDeadInstr();
    void    CalcDeadInstr(int step, int bot, int top, double * const * price);
    
    bool    UpIsOut; // this may be in a tweak ! so real instrument may not be dead
    bool    DownIsOut; // this may be in a tweak ! so real instrument may not be dead
    // below are barriers for each step
    CDoubleArray    UBarInput; // unadjsuted barrier levels
    CDoubleArray    LBarInput; // unadjsuted barrier levels
    CDoubleArray    UBarRebate; // rebate amount for upper barrier
    CDoubleArray    LBarRebate; // rebate amount for lower barrier
    double          UpperBar; // upper barrier at a tree slice after discrete monitor adjustment
    double          LowerBar; // lower barrier at a tree slice after discrete monitor adjustment
    // copied here for modification
    TBarrier        UType;
    TBarrier        LType;
    TBarDependence  BarDepend;
    TMonDependence  MonDepend;
    TPayoffMode     PayMode;

    int             NumOfPriceArray;

    bool            hasDefPayoff;   // used for dde

    // num of vol days per day, used for once a day barrier adjustment
    double          VolDaysPerYear;
    
    vector<double>  stepStrike;
    
private:
    const CDblBarrier* inst;
    vector<bool>       stepCanExercise;
    vector<bool>       stepCanKOUp;
    vector<bool>       stepCanKODown;
    bool               isStaticSpread;

    void SmoothPrice(const double* s, int step, int bot, int top, int pStart, int pEnd, const vector< double * > & p);
};

//////////// ***************************** /////////////////////////////////////////////
//     private class for all tree/FD product
//     new state variable interface for any num of factors
/////////////////////////////////////////////////////////
class DblBarrierFDProd : public LatticeProdEDRIns{
protected:
    typedef enum{CALL, PUT, BINARY, FORWARD} TPayoffMode;
    typedef enum{NA=-1, KI=0, KO=1} TBarrier;
    typedef enum{KI_KEEP_KO, KI_CANCEL_KO, ONCE_TOUCH, TWO_TOUCH} TBarDependence;
    typedef enum{BOTH, UPPER, LOWER} TMonDependence;

    const CDblBarrier*     inst;
    vector<double>      stepStrike;
    vector<bool>        stepCanExercise;
    vector<bool>       stepCanKOUp;
    vector<bool>       stepCanKODown;

    bool UpIsOut; // this may be in a tweak ! so real instrument may not be dead
    bool DownIsOut; // this may be in a tweak ! so real instrument may not be dead
    bool isDead;
    // below are barriers for each step
    CDoubleArray    UBarInput; // unadjsuted barrier levels
    CDoubleArray    LBarInput; // unadjsuted barrier levels
    CDoubleArray    UBarRebate; // rebate amount for upper barrier
    CDoubleArray    LBarRebate; // rebate amount for lower barrier
    double          UpperBar; // upper barrier at a tree slice after discrete monitor adjustment
    double          LowerBar; // lower barrier at a tree slice after discrete monitor adjustment
    // copied here for modification
    TBarrier        UType;
    TBarrier        LType;
    TBarDependence  BarDepend;
    TMonDependence  MonDepend;
    TPayoffMode     PayMode;

    // to allow switching between original and using slice operators update
    typedef void ( DblBarrierFDProd::*prod_FUNC )(
        int step, const TreeSlice & s, const vector< TreeSliceSP > & price );
    prod_FUNC prod_BWD_T, prod_BWD;

public:
    DblBarrierFDProd(const CDblBarrier* inst, FDModel* m);
    ~DblBarrierFDProd() {}

    // !!! this is to work with current tree which uses this to set up timeline only */
    virtual CVolRequestConstSP GetLNRequest() const;

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
    {
        return inst->ccyTreatment;
    }

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    /** initialisation, called ONCE only before initModel() for each new model instance */
    virtual void init(Control*  control) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing 
    virtual void initProd();

    /** calculate barriers and place barriers at inserted node if needed */
    virtual void preCalc(int step);

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type);

    /** check if knocked out already */
    void prepDeadInstr();
    /** calculate dead value */
    void calcDeadInstr( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price );

    // local methods
    /** product payoff method at maturity */
    void prod_BWD_T_orig( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price );
    void prod_BWD_T_oper( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price );

    // this is payoff boundary condition, for KO, early exercise etc.
    void prod_BWD_orig( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price );
    void prod_BWD_oper( int step, const TreeSlice & spot, const vector< TreeSliceSP > & price );

    // scale for notional 
    double scalePremium(const double& fairValue, YieldCurveConstSP disc);

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results);
};

DRLIB_END_NAMESPACE
#endif
