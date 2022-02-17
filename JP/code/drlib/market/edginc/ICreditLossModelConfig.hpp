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

#ifndef QLIB_ICREDITLOSSMODELCONFIG_HPP
#define QLIB_ICREDITLOSSMODELCONFIG_HPP

#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE_REF_COUNT(ICreditLossGen);
FORWARD_DECLARE(ICreditLossConfig);
FORWARD_DECLARE(IModelConfigMapper);
FORWARD_DECLARE(IModelConfigMapper);
FORWARD_DECLARE(CounterPartyCredit);

class MARKET_DLL ICreditLossModelConfig : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    ICreditLossModelConfig();
    virtual ~ICreditLossModelConfig();
    
    /**
     * Create a ICreditLossGen for the supplied ICreditLossConfig
     * object. The object returned might, for example, depend upon
     * the name of the ICreditLossConfig. Implementations are
     * expected to define an interface that the ICreditLossGen must
     * implement in order for this ICreditLossModelConfig to talk to
     * the ICreditLossConfig. Callers of this method are expected to know
     * the type of ICreditLossGen they are expecting to get back and
     * to dynamically cast to that type.
     * 
     * The "mapper" defines a way to retrieve which ICreditLossModelConfig
     * "inner models" to use in the case of recursive ICreditLossConfig.
     * "mapper" can be NULL, in which case it is up to this ICreditLossModelConfig
     * to define a way to retrieve ICreditLossModelConfig "inner models".
     * 
     * Note: there is no explicit mechanism to specify the expected return type of 
     * this method. In some cases, it means that the returned object will have to be
     * "poly-interfaces" to support different dynamic cast.
     * Example: the "lossGenerator" method applied on a tranche loss config should
     * return a "ICreditLossGenSP" that could be cast either to a "IEffectiveLossCurveGen"
     * (when pricing a simple tranche index) or a "ILossDistributionsGen" (when pricing a
     * CDO2).
     * */
    virtual ICreditLossGenSP lossGenerator(
        ICreditLossConfigConstSP lossConfig,
        IModelConfigMapperConstSP mapper) const = 0;

    virtual ICreditLossGenSP effCurveGenerator(
        ICreditLossConfigConstSP  lossConfig,
        CounterPartyCreditConstSP cpty,
        const bool                recoverNotional,
        IModelConfigMapperConstSP mapper) const = 0;

private:    
    ICreditLossModelConfig(const ICreditLossModelConfig& rhs); // don't use
    ICreditLossModelConfig& operator=(const ICreditLossModelConfig& rhs); // don't use
};

DECLARE(ICreditLossModelConfig);

DRLIB_END_NAMESPACE

#endif
