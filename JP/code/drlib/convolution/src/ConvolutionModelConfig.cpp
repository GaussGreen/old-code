#include "edginc/config.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/ConvolutionModelConfig.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/TrancheConvolutionLossGen.hpp"
#include "edginc/PortfolioNameConvolutionLossGen.hpp"
#include "edginc/NToDefaultLossConfig.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/FtDEffectiveCurveGen.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(ConvolutionModelConfigArray);

/** TYPE (for reflection) */        
CClassConstSP const ConvolutionModelConfig::TYPE =
CClass::registerClassLoadMethod("ConvolutionModelConfig",
                                typeid(ConvolutionModelConfig),
                                ConvolutionModelConfig::load);

/** Virtual destructor */
ConvolutionModelConfig::~ConvolutionModelConfig() {}

/** Constructor */
ConvolutionModelConfig::ConvolutionModelConfig():CObject(TYPE) {}

/** Explicit constructor */
ConvolutionModelConfig::ConvolutionModelConfig(
    IConvolutorSP convolutor,
    IConditionalDefaultsModelSP condDefaultsModel,
    bool returnBinaryDistribution):
        CObject(TYPE),
        convolutor(convolutor),
        returnBinaryDistribution(returnBinaryDistribution),
        condDefaultsModel(condDefaultsModel){}

IObject* ConvolutionModelConfig::defaultConstructor()
{
    return new ConvolutionModelConfig();
}

void ConvolutionModelConfig::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(ConvolutionModelConfig, clazz);
    IMPLEMENTS(ICreditLossModelConfig);
    IMPLEMENTS(IModelConfigMapper);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(convolutor, "Convolution engine");
    FIELD(returnBinaryDistribution,
        "Flag to control if we want to approximate the returned"
        "loss distribution by a binary distribution");
    FIELD(condDefaultsModel, "Conditional defaults model");
}

/** Dummy implementation of IModelConfigMapper (just returns "this") */
ICreditLossModelConfigConstSP ConvolutionModelConfig::innerModel(
    ICreditLossConfigConstSP lossConfig) const
{
    return ICreditLossModelConfigConstSP(this);
}

ICreditLossGenSP ConvolutionModelConfig::lossGenerator(
    ICreditLossConfigConstSP lossConfig,
    IModelConfigMapperConstSP mapper) const
{
    // TODO - use some double dispatch mechanism here
    if (CreditTrancheLossConfig::TYPE->isInstance(lossConfig))
    {
        const CreditTrancheLossConfig* tranche =
            DYNAMIC_CONST_CAST(CreditTrancheLossConfig, lossConfig.get());

        return ICreditLossGenSP(new TrancheConvolutionLossGen(
            tranche,
            tranche->trancheOutstandingNotional(false),
            tranche->getPortfolio().get(),
            mapper.get(),
            convolutor.get(),
            returnBinaryDistribution,
            condDefaultsModel.get(),
            tranche->getToday()));
	} 

	else if (FlatCDO2LossConfig::TYPE->isInstance(lossConfig))
    {
        const FlatCDO2LossConfig* slices =
			dynamic_cast<const FlatCDO2LossConfig*>(lossConfig.get());
            // DYNAMIC_CONST_CAST(FlatCDO2LossConfig, lossConfig.get());

		slices->getPortfolio().get();

		slices->trancheOutstandingNotional(false);

		condDefaultsModel.get();

		mapper.get();

		convolutor.get();

        return ICreditLossGenSP(new TrancheConvolutionLossGen(
            slices,
            slices->trancheOutstandingNotional(false),
            slices->getPortfolio().get(),
            mapper.get(),
            convolutor.get(),
            returnBinaryDistribution,
            condDefaultsModel.get(),
			slices->getToday()
		));
    }

    else if (PortfolioName::TYPE->isInstance(lossConfig))
    {
        return ICreditLossGenSP(new PortfolioNameConvolutionLossGen(
            DYNAMIC_CONST_CAST(PortfolioName, lossConfig.get()),
            condDefaultsModel.get()));
    }
    else
    {
        throw ModelException(
            "ConvolutionModelConfig does not handle " +
            lossConfig->getClass()->getName());
    }
}


/** Creates an IEffectiveCurveGen for an 'NtD' */
ICreditLossGenSP ConvolutionModelConfig::effCurveGenerator(
    ICreditLossConfigConstSP  lossConfig,
    CounterPartyCreditConstSP cpty,
    const bool                recoverNotional,
    IModelConfigMapperConstSP modelConfigMapper) const
{
    static const string method("ConvolutionModelConfig::effCurveGenerator");
    
    if (NToDefaultLossConfig::TYPE->isInstance(lossConfig)) {
        NToDefaultLossConfigConstSP ntdLossCfg(
            DYNAMIC_CONST_CAST(NToDefaultLossConfig, lossConfig.get()));

        return ICreditLossGenSP(new FtDEffectiveCurveGen(
            ntdLossCfg,
            modelConfigMapper,
            condDefaultsModel));
    }
    else {
        throw ModelException(method, "ConvolutionModelConfig does not handle " +
            lossConfig->getClass()->getName());
    }
}


/* external symbol to allow class to be forced to be linked in */
bool ConvolutionModelConfigLoad(){
    return (ConvolutionModelConfig::TYPE != 0);
}

DRLIB_END_NAMESPACE
