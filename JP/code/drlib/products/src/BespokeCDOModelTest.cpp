#include "edginc/config.hpp"
#include "edginc/BespokeCDOModelTest.hpp"
#include "edginc/FixedTrancheLossCalculator.hpp"
#include "edginc/ILossDistributionsGen.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/FunctionIntegrator1D.hpp"
#include "edginc/RectangleIntegrator1D.hpp"
#include "edginc/BetaConvolutor.hpp"
#include "edginc/RecursiveConvolutor.hpp"
#include "edginc/GaussianMarketFactorModel.hpp"
#include "edginc/CreditMetricsDefaultsModel.hpp"
#include "edginc/NToDefaultLossConfig.hpp"
#include "edginc/ConvolutionModelConfig.hpp"

DRLIB_BEGIN_NAMESPACE

// Utility (hacky) method to create the mapper "on the fly"
// according to "legacy" model parameters
void BespokeCDOModelTest::buildMapper(const bool useFastConvolution) const 
{
    // by defaults, build a "Credit Metrics" modelConfigMapper

    IIntegratorSP integrator(new FunctionIntegrator1D(
        Integrator1DSP(new RectangleIntegrator1D(101))));
    
    IConvolutorSP convolutor;
    
    if (useFastConvolution) {
        convolutor.reset(new BetaConvolutor(  
            1e-10,     // maximum variance value in % to be considered to be 0
            100,       // maximum nb of iterations in the beta function algorithm
            3.0e-7,    // precision in the beta function algorithm
            1.0e-30)); // flooring value of the internal variables close to 0 
    }
    else {
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

    modelConfigMapper.reset(new ConvolutionModelConfig(
        convolutor,
        condDefaultsModel,
        returnBinaryDistribution));
}

BespokeCDOModelTest::BespokeCDOModelTest(const CClassConstSP& clazz):
    CreditMetricsModel(clazz), modelConfigMapper(0){}

/** Destructor */
BespokeCDOModelTest::~BespokeCDOModelTest() {
}

class BespokeCDOModelTestHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BespokeCDOModelTest, clazz);
        SUPERCLASS(CreditMetricsModel);
        EMPTY_SHELL_METHOD(defaultBespokeCDOModelTest);
        FIELD(modelConfigMapper,
            "Mapper defining which model config to use for "
            "a given loss config.");
        FIELD_MAKE_OPTIONAL(modelConfigMapper);
    }
    
    static IObject* defaultBespokeCDOModelTest(){
        return new BespokeCDOModelTest(BespokeCDOModelTest::TYPE);
    }
};

CClassConstSP const BespokeCDOModelTest::TYPE = CClass::registerClassLoadMethod(
    "BespokeCDOModelTest", typeid(BespokeCDOModelTest), BespokeCDOModelTestHelper::load);

bool BespokeCDOModelTestLoad() {
    return (BespokeCDOModelTest::TYPE != 0);
}

class BespokeCDOModelTest::BespokeCDOModelLossCalculator:
    public virtual IFixedTrancheLossCalculator
{
public:
    
    // Constructor
    BespokeCDOModelLossCalculator(
        const DateTimeArray& timeline,
        IModelConfigMapperConstSP modelConfigMapper,
        ICreditLossConfigConstSP lossConfig):
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

void BespokeCDOModelTest::createLossCalculators(
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
    // Just in case it is not provided, init the modelConfigMapper
    initModelConfigMapper(cpty, recoverNotional, useFastConvolution(tranche));
    
    lossCalculator.reset(
        new BespokeCDOModelLossCalculator(timeline, modelConfigMapper, tranche));
}


/** Creates an IEffectiveCurveGen for an 'NtD' */
ICreditLossGenSP BespokeCDOModelTest::createEffCurveGenerator(
    NToDefaultLossConfigConstSP ntdLossCfg,
    CounterPartyCreditConstSP   cpty,
    const bool                  recoverNotional) const 
{
    // Just in case it is not provided, init the modelConfigMapper
    initModelConfigMapper(cpty, recoverNotional);

    // Retrieve model config
    ICreditLossModelConfigConstSP lossModelConfig =
        modelConfigMapper->innerModel(ntdLossCfg);

    return lossModelConfig->effCurveGenerator(ntdLossCfg,
                                              cpty,
                                              recoverNotional,
                                              modelConfigMapper);
}


// Builds a default for the modelConfigMapper if its current value is NULL
void BespokeCDOModelTest::initModelConfigMapper(
    CounterPartyCreditConstSP cpty,
    bool                      recoverNotional,
    const bool                useFastConvolution) const 
{
    static const string method("BespokeCDOModelTest::initModelConfigMapper");
    
    if (recoverNotional) {
        throw ModelException(method,
            "'Recover Notional' not yet supported.");
    }
    if (cpty.get() != 0) {
        throw ModelException(method,
            "Counterparty Credit Charge calculation not yet supported.");
    }

    if (!modelConfigMapper) {
        buildMapper(useFastConvolution);
    }
}

void BespokeCDOModelTest::createLossCalculatorsFIX(
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
		FlatCDO2LossConfigConstSP           ts) const           /* (I) */
{
	if (modelConfigMapper.get() == 0)
	{
		buildMapper(false);
	}
	
	lossCalculator.reset(
		new BespokeCDOModelLossCalculator(
			timeline,
			modelConfigMapper,
			ts)); 
};

DRLIB_END_NAMESPACE



