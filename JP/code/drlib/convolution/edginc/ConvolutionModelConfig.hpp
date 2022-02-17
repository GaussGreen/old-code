//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 04-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CONVOLUTIONMODELCONFIG_HPP
#define QLIB_CONVOLUTIONMODELCONFIG_HPP

#include "edginc/IModelConfigMapper.hpp"
#include "edginc/IConvolutor.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Model config implementation for loss configs from which we
 * can get a loss distribution using a convolution algorithm
 * */
class CONVOLUTION_DLL ConvolutionModelConfig:
    public CObject,
    public virtual ICreditLossModelConfig,
    public virtual IModelConfigMapper
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~ConvolutionModelConfig();

    /** Implementation of ICreditLossModelConfig */
    virtual ICreditLossGenSP lossGenerator(
        ICreditLossConfigConstSP lossConfig,
        IModelConfigMapperConstSP mapper) const;

    /** Creates an IEffectiveCurveGen for an 'NtD' */
    virtual ICreditLossGenSP effCurveGenerator(
        ICreditLossConfigConstSP  ntdLossCfg,
        CounterPartyCreditConstSP cpty,
        const bool                recoverNotional,
        IModelConfigMapperConstSP modelConfigMapper) const;

    /** Dummy implementation of IModelConfigMapper (just returns "this") */
    virtual ICreditLossModelConfigConstSP innerModel(
        ICreditLossConfigConstSP lossConfig) const;

    /** Explicit constructor */
    ConvolutionModelConfig(
        IConvolutorSP convolutor,
        IConditionalDefaultsModelSP condDefaultsModel,
        bool returnBinaryDistribution);

    // Return conditional defaults model
    IConditionalDefaultsModelSP CondDefaultsModel() { return condDefaultsModel; }
       
private:
    ConvolutionModelConfig(const ConvolutionModelConfig& rhs); // copy constructor - not defined
    ConvolutionModelConfig& operator=(const ConvolutionModelConfig& rhs); // operator= method - not defined

    ConvolutionModelConfig();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    /** Convolution engine */
    IConvolutorSP convolutor;

    // Flag to control if we want to approximate the returned
    // loss distribution by a "binary distribution" (in which case
    // we don't need to do a full convolution, we can just "convolute and
    // integrate").
    bool returnBinaryDistribution;

    /** Conditional default model */
    IConditionalDefaultsModelSP condDefaultsModel;
};

DECLARE(ConvolutionModelConfig);

DRLIB_END_NAMESPACE

#endif /*CONVOLUTIONMODELCONFIG*/
