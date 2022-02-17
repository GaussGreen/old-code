//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 14-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef EDR_BESPOKECDOMODELTEST_HPP
#define EDR_BESPOKECDOMODELTEST_HPP

#include "edginc/CreditMetricsModel.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(NToDefaultLossConfig);

/**
 * Temporary hack to test models for Bespoke CDO within the new
 * "RFL framework"
 * */
class PRODUCTS_DLL BespokeCDOModelTest: public CreditMetricsModel {
public:
    friend class BespokeCDOModelTestHelper;
    
    static CClassConstSP const TYPE;
    
    virtual ~BespokeCDOModelTest();

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

    /** Creates an IEffectiveCurveGen for an 'NtD' */
    virtual ICreditLossGenSP createEffCurveGenerator(
        NToDefaultLossConfigConstSP ntdLossCfg,
        CounterPartyCreditConstSP   cpty,
        const bool                  recoverNotional) const;

    // Return model config mapper 
    IModelConfigMapperSP ModelConfigMapper() { return modelConfigMapper; }

protected:
    /** Simple constructor */
    BespokeCDOModelTest(const CClassConstSP& clazz);

private:
    /** Don't use copy constructor */
    BespokeCDOModelTest(const BespokeCDOModelTest &rhs);
    
    /** Don't use */
    BespokeCDOModelTest& operator=(const BespokeCDOModelTest& rhs);
    
    class BespokeCDOModelLossCalculator; // define in source file
    
    // Utility (hacky) method to create the mapper "on the fly"
    // according to "legacy" model parameters
    void buildMapper(const bool useFastConvolution) const;
    
    // Builds a default for the modelConfigMapper if its current value is NULL
    void initModelConfigMapper(CounterPartyCreditConstSP cpty,
                               const bool                recoverNotional,
                               const bool                fastConvolution = true) const;

    // mutable so we can set it in buildMapper
    mutable IModelConfigMapperSP modelConfigMapper;
};

DRLIB_END_NAMESPACE
#endif //EDR_BESPOKECDOMODELTEST_HPP
