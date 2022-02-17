#include "edginc/config.hpp"
#include "edginc/BespokeCDOModelBCTest.hpp"
#include "edginc/FixedTrancheLossCalculator.hpp"
#include "edginc/ILossDistributionsGen.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/FunctionIntegrator1D.hpp"
#include "edginc/RectangleIntegrator1D.hpp"
#include "edginc/BetaConvolutor.hpp"
#include "edginc/RecursiveConvolutor.hpp"
#include "edginc/GaussianMarketFactorModel.hpp"
#include "edginc/CreditMetricsDefaultsModel.hpp"
#include "edginc/BCConvolutionModelConfig.hpp"

DRLIB_BEGIN_NAMESPACE

BespokeCDOModelBCTest::BespokeCDOModelBCTest(const CClassConstSP& clazz):
    CreditMetricsBaseCorrelation(clazz), modelConfigMapper(0){}

// Utility (hacky) method to create the mapper "on the fly"
// according to "legacy" model parameters
void BespokeCDOModelBCTest::buildMapper(
    YieldCurveConstSP discount, bool useFastConvolution) const
{
    // by defaults, build a "Credit Metrics" modelConfigMapper

    IIntegratorSP integrator(new FunctionIntegrator1D(
        Integrator1DSP(new RectangleIntegrator1D(101))));
    
    IConvolutorSP convolutor;
    
    if (useFastConvolution)
    {
        convolutor.reset(new BetaConvolutor(  
            1e-10,     // maximum variance value in % to be considered to be 0
            100,       // maximum nb of iterations in the beta function algorithm
            3.0e-7,    // precision in the beta function algorithm
            1.0e-30)); // flooring value of the internal variables close to 0 
    }
    else
    {
        convolutor.reset(new RecursiveConvolutor(
            5000, // max nb of slices for loss distribution
            1,    // nb by which to divide the GCD of notionals
            -1,   // loss unit override
            3));   // discretisation type
    }
        
    IMarketFactorModelSP marketFactorModel(new GaussianMarketFactorModel());
    marketFactorModel->validatePop2Object();

    bool returnBinaryDistribution = true;

    IConditionalDefaultsModelSP condDefaultsModel(
        new CreditMetricsDefaultsModel(
            integrator,
            marketFactorModel));

    CDoubleSP strikeMappingOverrideSP;
    if (strikeMappingOverride != 0)
    {
        strikeMappingOverrideSP.reset(CDouble::create(*strikeMappingOverride));
    }

    CDoubleSP spreadRatioOverrideSP;
    if (spreadRatioOverride != 0)
    {
        spreadRatioOverrideSP.reset(CDouble::create(*spreadRatioOverride));
    }

    BCConvolutionModelConfigSP bcTrancheModel(new BCConvolutionModelConfig(
        condDefaultsModel,
        convolutor,
        returnBinaryDistribution,
        compressionRatio,
        strikeMappingOverrideSP,
        lowerBCBetasOverride,
        upperBCBetasOverride,
        bcBetasTimeline,
        spreadRatioOverrideSP,
        spreadRatiosRegionalOverride,
        spreadRatiosRegions,
        historicalBetaAdjustment,
        namesOverride,
        indexWeightsOverride,
        authoriseNegativeEL,
        useExpectedLossRatio));

    // Hack to set up the discount curve - TODO remove ?
    bcTrancheModel->setDiscountCurve(discount);

    modelConfigMapper = bcTrancheModel;
}

/** Destructor */
BespokeCDOModelBCTest::~BespokeCDOModelBCTest() {
}

class BespokeCDOModelBCTestHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BespokeCDOModelBCTest, clazz);
        SUPERCLASS(CreditMetricsBaseCorrelation);
        EMPTY_SHELL_METHOD(defaultBespokeCDOModelBCTest);
        FIELD(modelConfigMapper,
            "Mapper defining which model config to use for "
            "a given loss config.");
        FIELD_MAKE_OPTIONAL(modelConfigMapper);
    }
    
    static IObject* defaultBespokeCDOModelBCTest(){
        return new BespokeCDOModelBCTest(BespokeCDOModelBCTest::TYPE);
    }
};

CClassConstSP const BespokeCDOModelBCTest::TYPE = CClass::registerClassLoadMethod(
    "BespokeCDOModelBCTest", typeid(BespokeCDOModelBCTest), BespokeCDOModelBCTestHelper::load);

bool BespokeCDOModelBCTestLoad() {
    return (BespokeCDOModelBCTest::TYPE != 0);
}

class BespokeCDOModelBCTest::BespokeCDOModelLossCalculator:
    public virtual IFixedTrancheLossCalculator
{
public:
    
    // Constructor
    BespokeCDOModelLossCalculator(
        const DateTimeArray& timeline,
        IModelConfigMapperConstSP modelConfigMapper,
        ICreditLossConfigConstSP lossConfig,
        YieldCurveConstSP discount):
            timeline(timeline),
            modelConfigMapper(modelConfigMapper),
            lossConfig(lossConfig)
    {
        // Retrieves model config
        ICreditLossModelConfigConstSP lossModelConfig =
            modelConfigMapper->innerModel(lossConfig);
        
        // Retrieves loss generator - this is where all the action related
        // to the model takes place
        ICreditLossGenSP lossGen =
            lossModelConfig->lossGenerator(lossConfig, modelConfigMapper);        
        
        // Casts to the expected loss generator
        ILossDistributionsGen* lossDistGen =
            dynamic_cast<ILossDistributionsGen*>(lossGen.get());
            
        // ALL COMPUTATION HERE
        IDistribution1DArraySP lossDist =
            lossDistGen->createLossDistributions(timeline);
    
        // and just store the results !!!
        losses.reset(new DoubleArray(timeline.size()));
        for (int i = 0; i < timeline.size(); ++i) {
			(*losses)[i] = (*lossDist)[i]->expectation();
		}
    }

    /** Calculate the expected loss for specified timepoint. */
    virtual void loss(
        int     timePoint,          // (I) do the calculation for this timepoint
        double& loss,                   /* (O) tranche loss amt  */
        double& lossCond) const    /* (O) tranche loss amt cond on cpty surviving */
    {
        loss = (*losses)[timePoint];
        lossCond = 0.0; // dummy value
    }
    
    /** store calculator specific results */
    virtual void storeResults(Results* result, Control* control) const
    {
        // do nothing...
    }
    
private:
    const DateTimeArray& timeline;
    IModelConfigMapperConstSP modelConfigMapper;
    ICreditLossConfigConstSP lossConfig;
    
    DoubleArraySP losses;
};

void BespokeCDOModelBCTest::createLossCalculators(
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
    IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc) const //(O)
{
    static const string method("BespokeCDOModelBCTest::createLossCalculators");
    
    if (recoverNotional)
    {
        throw ModelException(method,
            "'Recover Notional' not yet supported.");
    }
    if (cpty.get() != 0)
    {
        throw ModelException(method,
            "Counterparty Credit Charge calculation not yet supported.");
    }
    
    if (modelConfigMapper.get() == 0)
    {
        buildMapper(discount, useFastConvolution(tranche));
    }
    
    lossCalculator.reset(
        new BespokeCDOModelLossCalculator(
            timeline,
            modelConfigMapper,
            tranche,
            discount));
};

void BespokeCDOModelBCTest::createLossCalculatorsFIX(
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
		FlatCDO2LossConfigConstSP         ts            /* (I) */) const
{
	if (modelConfigMapper.get() == 0)
	{
		buildMapper(discount, false);
	}
	
	lossCalculator.reset(
		new BespokeCDOModelLossCalculator(
			timeline,
			modelConfigMapper,
			ts,
			discount)); 
};


DRLIB_END_NAMESPACE



