//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BasisZeroBondSeries.hpp
//
//   Description : Product used for testing port of basis diffusion from Rates SRM3
//
//   Date        : October 2005
//
//----------------------------------------------------------------------------

#ifndef BASIS_ZERO_BOND_SERIES_HPP
#define BASIS_ZERO_BOND_SERIES_HPP

#include "edginc/config.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/MCBasisCoupon.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE

/** BasisZeroBondSeries product 
    TO DO: UPDATE DESCRIPTION:
    Series of risky bonds -- pay one if survival and R if default
    valuation at "bondStart" for bond from "bondStart" to "bondMat" 
    valuation conditional on survival until bondStart */
class PRODUCTS_DLL BasisZeroBondSeries : public CInstrument,
                       public virtual MonteCarlo::IIntoProduct {
    
public:
    static CClassConstSP const TYPE;

    virtual ~BasisZeroBondSeries(){}

    void validatePop2Object();

    void GetMarket(const IModel*        model, 
                   const CMarketDataSP  market);
    
    virtual void Validate(); // after GetMarket
    
    DateTime getValueDate() const;

    /** Returns the name of the instrument's discount currency. */
    string discountYieldCurveName() const;

    /** Implementation of MonteCarlo::IntoProduct interface - the implementation of this is below */
    virtual MCProduct* createProduct(const MonteCarlo* model) const; 
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

protected:
    BasisZeroBondSeries(CClassConstSP clazz):
         CInstrument(clazz) {} // extra constructor
        
private:

    BasisZeroBondSeries(): CInstrument(TYPE) {}
    BasisZeroBondSeries(const BasisZeroBondSeries& rhs);            // not implemented
    BasisZeroBondSeries& operator=(const BasisZeroBondSeries& rhs); // not implemented

    static IObject* defaultBasisZeroBondSeries();
 
    /* fields */
    DateTime                valueDate;
    
    // BSERIES_DATA
    DateTimeArray           resetDate;      
    DateTimeArray           paymentDate; 
    
    // SP_INPUT
    /* base ir information 
    long      BaseIR_id;
    long      LiborTenor;   // underlying libor tenor in months
    char      LiborDCC;     // underlying libor DCC             
    char      LiborFreq;    // underlying libor Freq            
    long      LiborCrv;     // curve id for libor      
    is replaced by */
    
    YieldCurveWrapper       liborCurve;           
    MaturityPeriodSP        liborTenor;  // eg. 6M
    string                  accrualDCC; // Possibly different than libor Crv
    YieldCurveWrapper       discountYieldCurve; // BaseIR_id          
    
    // Internal
    DayCountConventionConstSP dcc;

    CAssetWrapper           fxAsset;    // irrelevant - to fix

    friend class BasisZeroBondSeriesMC;
};

class PRODUCTS_DLL BasisZeroBondSeriesMC : public MCProductClient,
                              public IMCStatelessProductClient {
public:
    // Ask for state variables 
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;
    
    BasisZeroBondSeriesMC(
        const BasisZeroBondSeries* inst,
        SimSeriesSP simSeries,
        InstrumentSettlementSP instSettle );

    virtual void payoff(
        const IPathGenerator* pathGen,
        Prices& prices);

    // IMCStatelessProductClient
    virtual IHistoricalContextSP createHistoricalContext();

    virtual IHistoricalContextSP getInitialHistoricalContext();

    virtual DateTimeArray getPastDates();

    virtual vector<int> finalize( const DateTimeArray& simDates );

    virtual void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history,
        MCProduct::Prices& prices );

private:
    // Update state variables 
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);
    DateTimeArray getTimeLine(set<int> idx) const;
    DateTimeArray getAllResetDates() const;

    double getDF( int t) const;

    double getEDF(int t,int paymentIdx) const;

    double getBasisCoupon(int t,int paymentIdx) const;

    struct PRODUCTS_DLL CashflowInfoHelper {
        CashflowInfoHelper(int idx) : paymentIdx(idx) {};
        int paymentDateIdx() const {
            return paymentIdx;
        }
//        bool operator < (const CashflowInfo & rhs) {
//            return paymentIdx < rhs.paymentIdx;
//        }
        int paymentIdx;
        int edfIdx; // offset inside associated SV's DateTimeArray
    };

    struct PRODUCTS_DLL SV {
        // Gives discount factor from t0 to the reset date
        SVDiscFactorSP df;
        // edf gives expected discount factors to all future payments dates
        SVExpectedDiscFactorSP edf;
        // basisCoupon gives basis coupon to all future payments dates
        MCBasisCoupon::IStateVarSP basisCoupon;
    };

    struct PRODUCTS_DLL GenSV {
        // Gives discount factor from t0 to the reset date
        SVGenDiscFactorSP df;
        // edf gives expected discount factors to all future payments dates
        SVGenExpectedDiscFactorSP edf;
        // basisCoupon gives basis coupon to all future payments dates
        MCBasisCouponSP basisCoupon;
    };
    
    struct PRODUCTS_DLL CashflowInfo {
        vector<CashflowInfoHelper> helper;
        SV sv;
        DateTime resetDate;
        set<int> getPaymentIndicies() const;
    };

    GenSV initializeGenSV(CashflowInfo& info); // modifies input

    vector<CashflowInfo> zeros;
    vector<GenSV> genSV;
    DateTimeArray allPaymentDates;
    
    const BasisZeroBondSeries* inst;
    map<DateTime,int> dateToIndex;     // Correspondence betweens key dates and index in sv array
    int lastResetIdx;

};

/* Extra fields in SP_INPUT to be captured elsewhere
typedef struct
{
    // deterministic yield curve 
    T_CURVE   BasisZCurve;

    // model and vol
    double    beta;
    double    q;
    double    spotVol;          // in bps vol terms. same size as Spread.
    double    BackboneCoeff;    // Back Bone Coefficient                  

    // Others 
    long      BasisType;    // specified as spread or percentage

} SP_INPUT;
*/
DRLIB_END_NAMESPACE

#endif  // BASIS_ZERO_BOND_SERIES_HPP
