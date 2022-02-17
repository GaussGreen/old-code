//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Jay Wang
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#ifndef QR_CREDITMETRICSABSCDOLOSSCALCULATOR_HPP
#define QR_CREDITMETRICSABSCDOLOSSCALCULATOR_HPP

#include "edginc/ABSCDOConvolution.hpp"
#include "edginc/CreditMetricsLossCalculatorBase.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/SkewSurface.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE_WRAPPER(CreditTrancheLossConfig);
FORWARD_DECLARE_WRAPPER(CounterPartyCredit);
class ConvolutionProduct;

class ISupportProxy {
public:
    virtual void getTopLoss(int timePoint, 
                            double& loss, 
                            double& lossCond) const = 0;
};

/** Implementation of ITrancheLossCalculator using 'CreditMetrics' */
class CONVOLUTION_DLL CreditMetricsABSCDOLossCalculator:
    public CreditMetricsLossCalculatorBase,
    public ISupportProxy 
{
public:
    virtual ~CreditMetricsABSCDOLossCalculator() {};

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CreditMetricsABSCDOLossCalculator(
        const DateTimeArray&           timeline,               /* (I) */
        CreditTrancheLossConfigConstSP tranche,                /* (I) */
        double                         lossUnit,               /* (I) */
        CounterPartyCreditConstSP      cpty,                   /* (I) */
        bool                           useSaddlePoint,         /* (I) */
        int                            numSaddlePoints,        /* (I) */
        int                            lossMarketFactors,      /* (I) */
        int                            decretionMarketFactors, /* (I) */
        bool                           fastConvolution,        /* (I) */
        const DoubleArray&             lossBetaOverride,       /* (I) */
        const DoubleArray&             decBetaOverride,        /* (I) */
        const DoubleArray&             lossDecBetaOverride,    /* (I) */
        bool                           calcExpectedNotional);  /* (I) */

    /** Calculate the expected loss for specified timepoint and strikes. 
        no caching supported yet */
    void loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double  k1,               /* (I) lower strike      */
        double  k2,               /* (I) upper strike      */
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const;  /* (O) tranche loss amt cond on cpty
                                     surviving */
    
    /** calculate tranche top loss */
    void getTopLoss(
        int     timePoint,
        double& loss,
        double& lossCond) const;

    /** store calculator results */
    void storeResults(Results* result, Control* control) const;

    /** Returns a key used to optimise repeated calculations of losses
        at the same timepoint. Not supported yet */
    virtual IKey* lossKey(
        int timePoint) const; // (I) do the calculation for this timepoint

    /** Same as above but allows the betas used to be overridden. The
        length of the array must be the same as the number of names or
        be empty (=> no beta overrides). Not supported yet */
    virtual IKey* lossKey(
        int                timePoint,   /* (I) do the calculation for 
                                           this timepoint */
        const DoubleArray& betaOverride) const; // (I) using these betas
        
private:
    void createMarketScenarios(const int lossMarketFactors, const int decretionMarketFactors);
    void createConditionalCurves(const DateTimeArray& timeline);
    
    void initializeOutputResults(const DateTimeArray&           timeline,
                                 CreditTrancheLossConfigConstSP tranche);

    void populateBasketInfoCM(CreditTrancheLossConfigConstSP     tranche,   /* (I) */
                              const DoubleArray&                 lossBetaOverride,
                              const DoubleArray&                 decBetaOverride,
                              ABSCDOConvolution::NameParamArray& basketInfo); // (O)

    void setupCM(const DoubleArray&        survivalProb, /* name survival proba */
                 double                    counterPartyProb) const; /* name survival proba */
    
    ABSCDOConvolution::NameParamSP createCptyInfoCM(
        CounterPartyCreditConstSP counterParty);
    
    /// fields ////
    const DateTimeArray&                timeline;
    double                              pastLoss;
    ABSCDOConvolution::NameParamArray   basketInfo;
    ABSCDOConvolution::NameParamSP      cpty;
    DoubleArray                         ntl;
    mutable DoubleArray                 density;
    mutable DoubleArray                 densityCond;
    int                                 maxs;
    int                                 maxl;
    double                              lossUnit;
    bool                                useSaddlePoint;
    int                                 numSaddlePoints;
    int                                 lossMarketFactors;
    int                                 decretionMarketFactors;
    bool                                useFastConvolution;
    bool                                calcExpectedNotional;
    // used to correlate loss and decretion market
    double                              lossMarketWeight;
    double                              decretionMarketWeight;
    
    GaussQuadIntegMethod                decIntMethod;
    GaussQuadIntegMethod                lossIntMethod;

    // output results, all these are mutable because loss method is const
    // 
    double                              initialNotional;
    double                              portNotional; // should be same as initialNotional except when you have past default
    mutable CashFlowArray               expectedTopLoss;
    mutable CashFlowArray               expectedBottomLoss;
    mutable CashFlowArray               expectedNotional;

    // MarketScenario Matrix
    ABSCDOConvolution::MarketScenarioMatrix   scenarios;
};

typedef smartPtr<CreditMetricsABSCDOLossCalculator> CreditMetricsABSCDOLossCalculatorSP;

class CONVOLUTION_DLL ProxyCalculator: public ITrancheLossCalculator {
public:
    ProxyCalculator(ISupportProxy* parent)
        : parent(parent) {}
    virtual ~ProxyCalculator() {};

    /** Calculate the expected loss for specified timepoint and strikes. 
        no caching supported yet */
    void loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double  k1,               /* (I) lower strike      */
        double  k2,               /* (I) upper strike      */
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const;  /* (O) tranche loss amt cond on cpty
                                     surviving */

    /** store calculator results */
    void storeResults(Results* result, Control* control) const {}

    /** Returns a key used to optimise repeated calculations of losses
        at the same timepoint. Not supported yet */
    virtual IKey* lossKey(
        int timePoint) const; // (I) do the calculation for this timepoint

    /** Same as above but allows the betas used to be overridden. The
        length of the array must be the same as the number of names or
        be empty (=> no beta overrides). Not supported yet */
    virtual IKey* lossKey(
        int                timePoint,   /* (I) do the calculation for 
                                           this timepoint */
        const DoubleArray& betaOverride) const; // (I) using these betas

private:
    ISupportProxy* parent;
};

// this really is just a combination of two CreditMetrics style calculators
class CONVOLUTION_DLL ABSCDOBaseCorrelationLossCalculator: 
    public ITrancheLossCalculator,
    public ISupportProxy 
{
public:
    virtual ~ABSCDOBaseCorrelationLossCalculator() {};

    ABSCDOBaseCorrelationLossCalculator(
        const DateTimeArray&           timeline,            /* (I) */
        CreditTrancheLossConfigConstSP tranche,             /* (I) */
        double                         lossUnit,            /* (I) */
        CounterPartyCreditConstSP      counterParty,        /* (I) */
        const DateTime&                maturity,           /* (I) */
        bool                           useSaddlePoint,      /* (I) */
        int                            numSaddlePoints,     /* (I) */
        int                            lossMarketFactors,   /* (I) */
        int                            decretionMarketFactors, /* (I) */
        bool                           fastConvolution,      /* (I) */
		bool							useLossSkew,		/* (I) */
		bool							useDecSkew,			/* (I) */
		bool							useLossDecSkew,		/* (I) */
        bool                           authoriseNegativeEL);

    /** Calculate the expected loss for specified timepoint and strikes. 
        no caching supported yet */
    void loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double  k1,               /* (I) lower strike      */
        double  k2,               /* (I) upper strike      */
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const;  /* (O) tranche loss amt cond on cpty
                                     surviving */

    /** calculate tranche top loss */
    void getTopLoss(
        int     timePoint,
        double& loss,
        double& lossCond) const;

    /** store calculator results */
    void storeResults(Results* result, Control* control) const;

    /** Returns a key used to optimise repeated calculations of losses
        at the same timepoint. Not supported yet */
    virtual IKey* lossKey(
        int timePoint) const; // (I) do the calculation for this timepoint

    /** Same as above but allows the betas used to be overridden. The
        length of the array must be the same as the number of names or
        be empty (=> no beta overrides). Not supported yet */
    virtual IKey* lossKey(
        int                timePoint,   /* (I) do the calculation for 
                                           this timepoint */
        const DoubleArray& betaOverride) const; // (I) using these betas
    
private:
    const DateTimeArray&                timeline;
    double                              initialNotional;
    
    double                              lowStrike;
    double                              highStrike;

    CreditMetricsABSCDOLossCalculatorSP calculatorLow; // for low strike
    CreditMetricsABSCDOLossCalculatorSP calculatorHigh; // for high strike

    // output requests
    mutable CashFlowArray               expectedTopLoss;
    mutable CashFlowArray               expectedBottomLoss;
    mutable CashFlowArray               expectedNotional;

    DoubleArraySP                       lowBetas;
    DoubleArraySP                       highBetas;

    bool                                authoriseNegativeEL;
};


DRLIB_END_NAMESPACE
#endif
