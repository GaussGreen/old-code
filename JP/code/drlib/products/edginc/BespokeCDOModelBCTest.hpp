//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 22-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef EDR_BESPOKECDOMODELBCTEST_HPP
#define EDR_BESPOKECDOMODELBCTEST_HPP

#include "edginc/CreditMetricsBaseCorrelation.hpp"
#include "edginc/IModelConfigMapper.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Temporary hack to test BC models for Bespoke CDO within the new
 * "RFL framework"
 * */
class PRODUCTS_DLL BespokeCDOModelBCTest:
    public CreditMetricsBaseCorrelation {
public:
    friend class BespokeCDOModelBCTestHelper;
    
    static CClassConstSP const TYPE;
    
    virtual ~BespokeCDOModelBCTest();

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

    virtual void createLossCalculatorsFIX(
        const DateTimeArray&                timeline,           /* (I) */    
        CounterPartyCreditConstSP           cpty,               /* (I) */
        const DateTime&                     maturity,           /* (I) */
        Control*                            control,            /* (I) */
        Results*                            results,            /* (I) */
        bool                                recoverNotional,    /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        IFixedTrancheLossCalculatorConstSP& lossCalculator,                   /* (O) */
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,      /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,              /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc,
		FlatCDO2LossConfigConstSP           ts) const;          /* (I) */

protected:
    /** Simple constructor */
    BespokeCDOModelBCTest(const CClassConstSP& clazz);

private:
    /** Don't use copy constructor */
    BespokeCDOModelBCTest(const BespokeCDOModelBCTest &rhs);
    
    /** Don't use */
    BespokeCDOModelBCTest& operator=(const BespokeCDOModelBCTest& rhs);
    
    class BespokeCDOModelLossCalculator; // define in source file

    // Utility (hacky) method to create the mapper "on the fly"
    // according to "legacy" model parameters
    void buildMapper(YieldCurveConstSP discount, bool useFastConvolution) const;
    
    // mutable so we can set it in buildMapper
    mutable IModelConfigMapperSP modelConfigMapper;
};

DRLIB_END_NAMESPACE
#endif //EDR_BESPOKECDOMODELBCTEST_HPP
