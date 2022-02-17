//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EqLinkCashFlow.hpp
//
//   Description   Coupon stream linked to Equity generalized performance
//
//
//   $Log: EqLinkCashFlow.hpp,v $
//----------------------------------------------------------------------------

#ifndef EDR_EQ_LINK_CASH_FLOW_HPP
#define EDR_EQ_LINK_CASH_FLOW_HPP

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/IAggregate.hpp"
#include "edginc/KOStubRule.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/TreeSlice.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************************************
// Allow some flexibility in styles of KOLeg by providing interface and 
// opportunity to extend choice of implementers. 
class PRODUCTS_DLL IKOLeg {
public:
    virtual double getLevel(int simDateIdx) const = 0;

    // if hit, when pay. NOT settle adjusted
    virtual DateTime getPayDate(const DateTime& hitDate) const = 0;

    virtual ~IKOLeg() {};
};

typedef refCountPtr<IKOLeg> IKOLegSP;

// Build IKOLeg via a Maker class.
// It's a bit like building the thing we need in 2 stages.
class PRODUCTS_DLL IKOLegMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    virtual IKOLegSP getKOLeg(const DateTimeArray& simDates,
                                const YieldCurve*    discount) = 0;

    // populate from market cache 
    virtual void getMarket(const IModel*     model, 
                           const MarketData* market) = 0;

    virtual ~IKOLegMaker() {};

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IKOLegMaker> IKOLegMakerSP;
*/

//////////////////////////////////////////////////////////////
// maker class
//////////////////////////////////////////////////////////////
class PRODUCTS_DLL     EqLinkCashFlow: public CObject {
public:
    static CClassConstSP const TYPE;  
    friend class EqCpnKOLeg;

    void validatePop2Object();
    
    EqLinkCashFlow(); 

    EqLinkCashFlow(const IDoubleArrayModifierMakerSP eqLinkPerf,
                    const DateTimeArray observDates);

    EqLinkCashFlow(const IDoubleArrayModifierMakerSP eqLinkPerf,
                   const DateTimeArray observDates,
                   const DateTimeArray paymentDates);

    DateTimeArray getObservDates() const;
    
	DateTimeArray  getPaymentDates() const;

    IDoubleArrayModifierMakerSP  getEqLinkPerf();

    CashFlowArray getKnownCashFlow(CashFlowArray samples, 
                                    DateTime valDate, 
                                    double refLevel,
                                    const InstrumentSettlementSP  instSettle,
                                    const CAssetWrapper    asset,
                                    bool isExcludeKnownToday = false);

private:
    /// Invoked when Class is 'loaded' 
    static void load(CClassSP& clazz);
    
    static IObject* defaultEqLinkCashFlow();

    // data
    IDoubleArrayModifierMakerSP  eqLinkPerf;        // performace of each observation date
    DateTimeArray                observDates;       // maturities of each eq link coupons.
    DateTimeArray                paymentDates;
};

typedef smartPtr<EqLinkCashFlow> EqLinkCashFlowSP;

//////////////////////////////////////////////////////////////
// maker class with KO rule
//////////////////////////////////////////////////////////////
class PRODUCTS_DLL     EqCpnKOMaker: public CObject {
public:
    static CClassConstSP const TYPE;  
    friend class EqCpnKOLeg;

    void validatePop2Object();
    
    EqCpnKOMaker(); 

    EqCpnKOMaker(const EqLinkCashFlowSP eqCpn,
                    const DateTime  initialAccrueDate,
                    const KOStubRuleSP koRule,
                    bool inAdvance = false);

    DateTimeArray getObservDates();
    
	DateTimeArray  getPaymentDates();

    IDoubleArrayModifierMakerSP  getEqLinkPerf();

    KOStubRuleSP getKOStubRule();

    DateTime      getInitialAccrueDate();

    CashFlowArray getKnownCashFlow(CashFlowArray samples, 
                                    DateTime valDate, 
                                    double refLevel,
                                    const InstrumentSettlementSP  instSettle,
                                    const CAssetWrapper    asset,
                                    bool isExcludeKnownToday = false);

    // return next payment date corresponding to hit Date.
    // currently, "S" is not supported.  So always payment date is next pay date.
    bool getPayDate(const DateTime hitDate,const InstrumentSettlementSP  instSettle,
                    DateTime* payDate) const;

private:
    /// Invoked when Class is 'loaded' 
    static void load(CClassSP& clazz);
    
    static IObject* defaultEqCpnKOMaker();

    // data
    EqLinkCashFlowSP eqCpn;
    KOStubRuleSP     koRule;
    DateTime         initialAccrueDate; // first accrue start Date.
    bool             inAdvance;
};

typedef smartPtr<EqCpnKOMaker> EqCpnKOMakerSP;

//////////////////////////////////////////////////////////////
//      Internal Class
//      EqKOLeg (with KO rule)
//////////////////////////////////////////////////////////////
class EqCpnKOLeg;
typedef refCountPtr<EqCpnKOLeg> EqCpnKOLegSP;

class PRODUCTS_DLL EqCpnKOLeg{
public:
    EqCpnKOLeg();
    
    EqCpnKOLeg(const EqLinkCashFlowSP els,
             SimpleDoubleArray& perfs,
             DateTime valueDate,
             KOStubRuleSP koRule,
             bool inAdvance,
             DateTime initialAccrueDate = DateTime(0,0));

    // constructor with EqCpnKOMaker
	EqCpnKOLeg(const EqCpnKOMakerSP elcKO,
                SimpleDoubleArray& perfs,
                DateTime valueDate);

    // This is for Tree.  In Tree, it's not easy to handle the coupon whose observation in ealier than
    // KO date, because it's necessary to manage with state varaible tree (time consuming).
    // Thus, for time being, "S" is not available.   
    // Now, it just return array of 1 (B) or 0 (N).
    // When the time step is observation date, then "B" : Do not cancel the cpn /  "N" : Cancel the cpn.
    // It also look at "inAdvance".  In case of inAdvance & B, it return 0 for observe date to be cancelled.
    // Not Only exact time point, the same date is also checked for Tree usage. 
    // (Need to moidify for MC use).
    void setKOFactor(const DateTimeArray timeLine, 
                     const vector<bool>& stepIsKO, 
                     vector<double>& factors);

    // set up the discount factor.  from payment date to observ date (mainly for tree).
    // check whether there is paydates or not, and judge use instSettle or not.
	void setDiscFacts(const DateTimeArray timeLine,
                      const YieldCurveConstSP disoucnt,
                      const InstrumentSettlementSP  instSettle,
                      const CAssetWrapper asset);

    // set flag whether the date step is monitoring date or not.
    IntArray setMonStep(const DateTimeArray timeLine);

	// set the payment date corresponding to monitoring dates.
    void setPayDates(const DateTimeArray timeLine,const InstrumentSettlementSP  instSettle);

    // set up the payment by using EqCashFlow's paymentDates
    void setPayDatesAndDiscFacts(const DateTimeArray timeLine, 
                                const YieldCurve* discount,
                                const DateTime toDate);

    // return the payment date corresponding to simulation dates.  
    // i.e. The size of array  is same to timeLine, and could have duplicate dates.
    DateTimeArray getPayDatesOnTimeLine() const ;

    // return the payment date of elc, after removing duplicate of input.
    DateTimeArray getPayDates() const;
    
    // return index of payment dates array (after removed duplicated one),
    // which corresponding to the timeLine[step]
    int           payDateIndex(int step);

    // return payment date,
    // which corresponding to the timeLine[step]
    DateTime payDateOnKO(int step);

    // return with PV adjust.
    double getValue(int step);

    class getValue_oper : public SliceMarker< getValue_oper >
    {
    public:
        getValue_oper(const TreeSlice & spot,
                      int step,
                      double refLevel,
                      double &eqPerf,
                      EqCpnKOLeg & eqCpn);

        // TreeSlice "expression template" primitives
        static const int sliceCount = 1;
        template< typename S >
        const S** listSlices(const S** list) const
        {
            return spot.listSlices( list ) ;
        }
        inline double calc() const
        {
            return apply( spot.calc() );
        }
        void printDebug(char *s) const
        {
            strcat(s, "(PenultSmooth)");
        }

    private:
        double apply( double s ) const;

        const TreeSlice & spot;
        int step;
        double refLevel;
        mutable double *pPerf;
        EqCpnKOLeg & eqCpn;
    };

    
    // for MC product, calculate the range coupon value.
	double getValue(const int   startStep,
                    const int   endStep,  
                    const int   nbAssets);

    // for MC SV product, calculate the range coupon value.  Use SV DiscFactors.
    double getValue(const int   startStep,            //(I) start step of path
                    const int   endStep,              //(I) calculate up to endStep                             
                    const int   nbAssets,
                    SVDiscFactorSP discFactsSV);

    void  makeKnownCashFlow(const int   startStep,
                            const int   endStep,  
                            const int   nbAssets,
                            const bool  hasFuture);

    // return a copy of knownCFL
    CashFlowArraySP getKnownCashFlows();

    // a class to get the genereric performance results.
    class PRODUCTS_DLL PerfContainer{
    public:
        PerfContainer();
        PerfContainer(IDoubleArray* components);

        double getPerf() ;

        double getPerf(int index) ;

    private:
        IDoubleArray*          components;
    };

public:
    bool        isNull;     // whether the data is initialized or not.
    EqLinkCashFlowSP elc;
    KOStubRuleSP  koRule;        

private:
    IDoubleArrayModifierSP    performance;
    PerfContainer       couponPayments;

    DateTime        initialAccrueDate;
//        IDoubleArrayModifierMakerSP    perfMaker;        
//        DateTimeArray observDates;

    // whether the observation (fixing) is inAdvance or not (inArrear).
    // this flag change the rule when the KO date is same to
    // observDates.  inAdvance case, the Cpn should be cancelled
    // but inAdvance = false case, the cpn shouldn't be cancelled.
    // it's similar to koStubRule, but this handle only the case 
    // when koDate is same to observe date.  i.e. < or <=.
    // koStubRule manage those cancelling rule more longer period.
    bool    inAdvance;

    DateTimeArray payDates;             // payment dates if the KO happens on simulation dates. (Size is same as time steps)
    DoubleArray   discFacts;            // store the PV factors.
    
    BoolArray         isPayStep;        // is this step determine the payoff of each coupn?
                                        // The array size is time steps of engine (MC simDates or Tree TimeSteps)
    IntArray          monStepMap;       // int array which tells you which index date is corresponding to eq observDates.
    IntArray          obsIndex;         // to know observeDate index which corresponding to a simulation time step.
    CashFlowArray     knownCFL;         // knownCFL

    IntArray          payIndex;         // map from observDates array to unique payment dates array

    DateTime          valueDate;        // Need to store to calculate not Paid Coupons.
    double            remainedValue;    // PV of CashFlow which are not settled.
    bool              isAlreadyApplied; // no need to performance->apply() anymore.

    // this is private function!!
    double getPerf(const int idx);

    // set up the discount factor.  from payment date to observ date (mainly for tree).
    // only available when payment dates are given.
	void setDiscFacts(const DateTimeArray timeLine,const YieldCurveConstSP discount);
    
    // for State Variable Usage....
    bool hasDiscFactsSV;
    SVDiscFactorSP discFactsSV;
};



DRLIB_END_NAMESPACE

#endif