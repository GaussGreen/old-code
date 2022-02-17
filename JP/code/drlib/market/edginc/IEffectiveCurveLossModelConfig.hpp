//----------------------------------------------------------------------------
//                                                                           
// Group       : CH Quantitative Research                                    
//                                                                           
// Description : Interface class to control which instance of a model is used
//               for a given ICreditLossConfig                               
//                                                                           
// Date        : July 2006                                                   
//                                                                           
//----------------------------------------------------------------------------

#ifndef QLIB_IEFFECTIVECURVELOSSMODELCONFIG_HPP
#define QLIB_IEFFECTIVECURVELOSSMODELCONFIG_HPP

#include "edginc/Object.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/ICreditLossModelConfig.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(ICreditLossConfig);
FORWARD_DECLARE(IEffectiveCurveLossModelConfig);
FORWARD_DECLARE(CreditTrancheLossConfig);
FORWARD_DECLARE(FlatCDO2LossConfig);
FORWARD_DECLARE_REF_COUNT(ICreditLossGen);
FORWARD_DECLARE_REF_COUNT(IEffectiveCurveLossGen);
FORWARD_DECLARE_REF_COUNT(IEffectiveCurveGen);
FORWARD_DECLARE(NToDefaultLossConfig);
FORWARD_DECLARE(CounterPartyCredit);
FORWARD_DECLARE(IModelConfigMapper);

class MARKET_DLL IEffectiveCurveLossModelConfig : 
    public virtual ICreditLossModelConfig 
{
public:
    static CClassConstSP const TYPE;

    IEffectiveCurveLossModelConfig();
    virtual ~IEffectiveCurveLossModelConfig();

    /** The interface that ICreditLossConfig objects must implement in order 
        for this EffectiveCurveLossModelConfig to be used to produce an
        effective curve loss generator */
    class MARKET_DLL IIntoLossGen {
    public:
        static CClassConstSP const TYPE;

        virtual ~IIntoLossGen();

        /** Create an IEffectiveCurveLossGen as specified by the supplied
         * IEffectiveCurveLossModelConfig */
        virtual IEffectiveCurveLossGenSP lossGenerator(
            IEffectiveCurveLossModelConfigConstSP effCurveLossModelConfig) const = 0;
    };

    /** The interface that ICreditLossConfig objects must implement in order 
        for this EffectiveCurveLossModelConfig to be used to produce an
        effective curve generator */
    class MARKET_DLL IIntoEffCurveGen {
    public:
        static CClassConstSP const TYPE;

        virtual ~IIntoEffCurveGen();

        /** Create an IEffectiveCurveGen as specified by the supplied
         * IEffectiveCurveLossModelConfig */
        virtual ICreditLossGenSP effCurveGenerator(
            IEffectiveCurveLossModelConfigConstSP effCurveLossModelConfig,
            CounterPartyCreditConstSP             cpty,
            const bool                            recoverNotional,
            IModelConfigMapperConstSP             mapper) const = 0;
    };

    /** Creates an IEffectiveCurveLossGen for a 'tranche'*/
    virtual IEffectiveCurveLossGenSP createLossGenerator(
        CreditTrancheLossConfigConstSP trancheLossCfg) const = 0;

    virtual IEffectiveCurveLossGenSP createLossGenerator(
        FlatCDO2LossConfigConstSP tsLossCfg) const = 0;

    /** Creates an IEffectiveCurveGen for an 'NtD' */
    virtual ICreditLossGenSP createEffCurveGenerator(
        NToDefaultLossConfigConstSP ntdLossCfg,
        CounterPartyCreditConstSP   cpty,
        const bool                  recoverNotional) const = 0;

private:
    IEffectiveCurveLossModelConfig(const IEffectiveCurveLossModelConfig& rhs); // don't use
    IEffectiveCurveLossModelConfig& operator=(const IEffectiveCurveLossModelConfig& rhs); // don't use
};


DRLIB_END_NAMESPACE

#endif
