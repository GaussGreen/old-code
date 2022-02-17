//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditMetricsABSCDO.hpp
//
//   Description : Credit Metrics applied to ABS CDO
//
//   Date        : Jan 2006
//
//----------------------------------------------------------------------------

#ifndef EDR_CREDIT_METRICS_ABSCDO_HPP
#define EDR_CREDIT_METRICS_ABSCDO_HPP

#include "edginc/CreditMetricsModel.hpp"
#include "edginc/SkewSurface.hpp"

DRLIB_BEGIN_NAMESPACE
class ITrancheLossCalculatorLegacy;
/** 
 * CM model for ABS CDO
 * */
class PRODUCTS_DLL CreditMetricsABSCDO : public CreditMetricsModel {
public:
    friend class CDO;

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~CreditMetricsABSCDO();
    
    virtual void createLossCalculators(
        const DateTimeArray&                timeline,           /* (I) */
        CreditTrancheLossConfigConstSP      tranche,            /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        const DateTime&                     maturity,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        bool                                recoverNotional,    /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        IFixedTrancheLossCalculatorConstSP& lossCalculator,               /* (O) */
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,  /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,          /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const; //(O)

    /** Invoked by the containing model before the instrument data is fetched
        ie before CInstrument::GetMarket is invoked. Overrides
        CreditMetricsModel in order to get correlation between 
        decretion and loss markets */
    virtual void preFetchMarketData(Model*            model,
                                    MarketDataConstSP market);

    /** Called immediately after object constructed */
    void validatePop2Object();


	virtual const bool hasStochasticRecoveries(
		CreditTrancheLossConfigConstSP tranche) const
	{ return true; }

	virtual const bool modelRecoveredNotional(
		CreditTrancheLossConfigConstSP tranche) const
	{ return true; }


protected:
    /** Only build instances of that class using reflection */
    CreditMetricsABSCDO(const CClassConstSP& clazz);
    
private:
    /** Default constructor */
    static IObject* defaultConstructor();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** For CreditMetrics::IIntoProduct */
    static void loadIntoProduct(CClassSP& clazz);

    bool useSaddlePoint;
    int numSaddlePoints;
    int lossMarketFactors;
    int decretionMarketFactors;

    bool useLossSkew;
	bool useDecSkew;
	bool useLossDecSkew;

    /**
     * Flag to authorise negative expected losses (=TRUE).
     * Otherwise negative expected losses will be floored to 0.
     * Default is FALSE (i.e. do NOT authorise negative expected losses).
     * */
    bool authoriseNegativeEL;

protected:    
    // ------
    // FIELDS
    // ------

};


DRLIB_END_NAMESPACE

#endif
